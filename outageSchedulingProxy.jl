#!~/Downloads/julia-0.4.5/usr/bin/julia

using ClusterManagers
using OffsetArrays
using ArgParse

include("./dataInputOutput.jl") #load this on the main process already
include("./STproxies.jl")

function outageSchedulingProxy(path, filereq, s0, s1, day0, day1, M, NbOfCPUs, NbOfNodes, version)

    if NbOfNodes > 0
        addprocs(SlurmManager(M), cpus_per_task=NbOfCPUs, nodes=NbOfNodes, time="2-00:00:00", mem_per_cpu=4000, exclude="node[120-128]")
    else
        addprocs(M)
    end

	#these includes come after addprocs() so files are loaded on all processes
	include("./dataInputOutput.jl")
	include("./STproxies.jl")

    filesch = "$(filereq[1:end-4])_Schedule_MS$(s0)-$(s1)_Days$(day0)-$(day1)_v$(version).csv"

    req = try readcsv(filereq, Int) catch; [] end
    if isempty(req) return end
    sch = try readcsv(filesch, Int) catch; [] end

    #NOTE:
    #req is a matrix where each line is a "pair" (branch_id, duration)
    #sch is a matrix where each line is a "triplet" (branch_id, starting_day, ending_day)

	while(true)

		sch_branches = sch[:,1]

		req_left = []
		for i=1:size(req,1)
			branch_id = req[i,1]
			if !in(branch_id, sch_branches)
				req_left = vcat(req_left, reshape(req[i,:], 1, 2))
			end
		end

		println("requests left = $(req_left[:,1])")
		if isempty(req_left) break end

		newsch = scheduleNext(path, req_left, sch, s0, s1, day0, day1, version)
        sch = vcat(sch, newsch)

        writecsv(filesch, sch)

	end

    for i in workers() rmprocs(i) end

end

function scheduleNext(path, req_left, sch, s0, s1, day0, day1, version)

	newsch = []
	maxCompoundImpact = typemin(Float64)
	for i=1:size(req_left,1)
		branch_id = req_left[i,1]
		duration = req_left[i,2]
		Impact = computeOutageImpact(path, s0, s1, day0, day1, branch_id, sch, version)
		AvgImpact = mean(Impact, 1)[1,:]
		CompoundImpact = OffsetArray(Float64, day0:day1 - duration + 1)
		for day=day0:day1 - duration + 1
			CompoundImpact[day] = sum(AvgImpact[day:day+duration-1])
		end
		if minimum(CompoundImpact) > maxCompoundImpact
			maxCompoundImpact = minimum(CompoundImpact)
            starting_day = day0 + indmin(CompoundImpact) - 1    #this is correct, since indmin() does not work as expected in OffsetArrays
            ending_day = starting_day + duration - 1
			newsch = [branch_id starting_day ending_day]
		end
	end

    return newsch

end

function computeOutageImpact(path, s0, s1, day0, day1, o, sch, version)

	outB = OffsetArray([[] for day=day0:day1], day0:day1)
    Impact = OffsetArray(Float64, s0:s1, day0:day1)

    #initialize outage schedule array outB:
	if !isempty(sch)
		for i=1:size(sch,1)
			branch_id = sch[i,1]
			starting_day = sch[i,2]
			ending_day = sch[i,3] 
			for day=starting_day:ending_day
				push!(outB[day], branch_id)
			end
		end
	end

    np = nprocs()  # determine the number of processes available
    s = s0
    day = day0
    # function to produce the next work item from the queue.
    # in this case it's just an index.
    nextidx() = (nexts=s; nextday=day; day+=1; if day > day1 day=day0; s+=1 end; (nexts, nextday))
    @sync begin
        for p=1:np
            if p != myid() || np == 1
                @async begin
                    while true
                        sidx, dayidx = nextidx()
                        if sidx > s1
                            break
                        end
                        println("$(sidx), $(dayidx)")
                        outBbase = outB[dayidx]
                        outBnew = copy(outBbase)
                        push!(outBnew, o)
                        if version == 1
                            BaseCost = 0
                        elseif version == 2
                            DAoutcome, RToutcome = remotecall_fetch(readOrComputeOutageOutcome, p, path, sidx, dayidx, outBbase)
                            BaseCost = DAoutcome[:DAcost][1] + sum(RToutcome[:RTcost])
                        end
                        DAoutcome, RToutcome = remotecall_fetch(readOrComputeOutageOutcome, p, path, sidx, dayidx, outBnew)
                        NewCost = DAoutcome[:DAcost][1] + sum(RToutcome[:RTcost])
                        Impact[sidx, dayidx] = NewCost - BaseCost
                    end
                end
            end
        end
    end   

	return Impact

end


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "path"
            help = "the path to the folder with the system specification files"
            arg_type = AbstractString
            default = "./data/IEEE-RTS96"
        "filereq"
            help = "the file with the outage requests"
            arg_type = AbstractString
            default = "./data/IEEE-RTS96/out_sch_cases/Case1OutageRequestData.csv"
        "M"
            help = "Number of tasks"
            arg_type = Int
            default = 4
        "NbOfCPUs"
            help = "Number of CPUs per task"
            arg_type = Int
            default = 0
        "NbOfNodes"
            help = "Number of nodes"
            arg_type = Int
            default = 0
        "s0"
            help = "The first micro-scenario"
            arg_type = Int
            default = 2
        "s1"
            help = "The last micro-scenario"
            arg_type = Int
            default = 3
        "day0"
            help = "The first day"
            arg_type = Int
            default = 1
        "day1"
            help = "The last day"
            arg_type = Int
            default = 56
        "version"
            help = "which version to use (1: absolute cost, or 2: relative cost)"
            arg_type = Int
            default = 1
    end

    return parse_args(s)
end

function main()

	parsed_args = parse_commandline()
	println("Parsed args:")
	for (arg,val) in parsed_args
	    println("  $arg  =>  $val")
	end

	path = parsed_args["path"]
	filereq = parsed_args["filereq"]
	M = parsed_args["M"]
	NbOfCPUs = parsed_args["NbOfCPUs"]
	NbOfNodes = parsed_args["NbOfNodes"]
	s0 = parsed_args["s0"]
	s1 = parsed_args["s1"]
	day0 = parsed_args["day0"]
	day1 = parsed_args["day1"]
    version = parsed_args["version"]


	outageSchedulingProxy(path, filereq, s0, s1, day0, day1, M, NbOfCPUs, NbOfNodes, version)

end

main()

