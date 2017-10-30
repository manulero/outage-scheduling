#!~/Downloads/julia-0.4.5/usr/bin/julia


using ClusterManagers
using OffsetArrays
using ArgParse

include("./dataInputOutput.jl") #load this on the main process already
include("./STproxies.jl")



function assessmentBaseCase(path, s0, s1, day0, day1, M, NbOfCPUs, NbOfNodes)

    if NbOfNodes > 0
	   addprocs(SlurmManager(M), cpus_per_task=NbOfCPUs, nodes=NbOfNodes, time="2-00:00:00", mem_per_cpu=4000, exclude="node[120-128]")
    else
        addprocs(M)
    end
    #these includes come after addprocs() so files are loaded on all processes
    include("./dataInputOutput.jl")
    include("./STproxies.jl")


    np = nprocs()  # determine the number of processes available
    s = s0
    day = day0
    # function to produce the next work item from the queue.
    # in this case it's just an index.
    #nextidx() = (nexts=s; nextday=day; day+=1; if day > day1 day=day0; s+=1 end; (nexts, nextday))
    nextidx() = (nexts=s; nextday=day; s+=1; if s > s1 s=s0; day+=1 end; (nexts, nextday))
    @sync begin
        for p=1:np
            if p != myid() || np == 1
                @async begin
                    while true
                        sidx, dayidx = nextidx()
                        if dayidx > day1
                            break
                        end
                        println("$(sidx), $(dayidx)")
                        BaseOutcome = remotecall_fetch(readOrComputeOutageAssessmentOutcome, p, path, sidx, dayidx, [])
                        #println("$(sidx), $(dayidx), $(BaseOutcome)")
                    end
                end
            end
        end
    end


    for i in workers() rmprocs(i) end

end


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "path"
            help = "the path to the folder with the system specification files"
            arg_type = AbstractString
            default = "./data/IEEE-RTS96"
        "M"
            help = "Number of tasks"
            arg_type = Int
            default = 0
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
            default = 1
        "s1"
            help = "The last micro-scenario"
            arg_type = Int
            default = 1
        "day0"
            help = "The first day"
            arg_type = Int
            default = 7
        "day1"
            help = "The last day"
            arg_type = Int
            default = 7
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
	M = parsed_args["M"]
	NbOfCPUs = parsed_args["NbOfCPUs"]
	NbOfNodes = parsed_args["NbOfNodes"]
	s0 = parsed_args["s0"]
	s1 = parsed_args["s1"]
	day0 = parsed_args["day0"]
	day1 = parsed_args["day1"]


	assessmentBaseCase(path, s0, s1, day0, day1, M, NbOfCPUs, NbOfNodes)

end

main()
