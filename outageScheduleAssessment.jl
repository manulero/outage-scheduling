
using OffsetArrays
using ArgParse

include("./dataInputOutput.jl") #load this on the main process already
include("./STproxies.jl")



function outageScheduleAssessment(path, filesch, s0, s1, day0, day1, M)

    addprocs(M)

    #these includes come after addprocs() so files are loaded on all processes
    include("./dataInputOutput.jl")
    include("./STproxies.jl")


    sch = try
        readcsv(filesch, Int)
    catch err
        println(err)
        []
    end

    if isempty(sch)
        return
    end


    outB = OffsetArray([[] for day=day0:day1], day0:day1)

    for i=1:size(sch,1)
        A = collect(day0:day1)
        B = collect(sch[i,2]:sch[i,3])
        C = intersect(A,B)
        if isempty(C) continue end
        branch_id = sch[i,1]
        starting_day = C[1]
        ending_day = C[end] 
        for day=starting_day:ending_day
            push!(outB[day], branch_id)
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
                        if isempty(outB[dayidx]) continue end
                        BaseCost = remotecall_fetch(readOrComputeOutageAssessment, p, path, sidx, dayidx, [])
                        NewCost = remotecall_fetch(readOrComputeOutageAssessment, p, path, sidx, dayidx, outB[dayidx])
                        println("$(sidx), $(dayidx), $(BaseCost), $(NewCost)")
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
        "filesch"
            help = "the file with the outage schedule"
            arg_type = AbstractString
            #default = "./data/IEEE-RTS96/out_sch_cases/Case1OutageRequestData_Schedule_MS1-10_Days1-182.csv"
            default = "./data/IEEE-RTS96/out_sch_cases/Case1OutageRequestData_Schedule_MS1-2_Days1-182_POWERTECH.csv"
        "M"
            help = "Number of tasks"
            arg_type = Int
            default = 4
        "s0"
            help = "The first micro-scenario"
            arg_type = Int
            default = 33
        "s1"
            help = "The last micro-scenario"
            arg_type = Int
            default = 48
        "day0"
            help = "The first day"
            arg_type = Int
            default = 1
        "day1"
            help = "The last day"
            arg_type = Int
            default = 182
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
	filesch = parsed_args["filesch"]
	M = parsed_args["M"]
	s0 = parsed_args["s0"]
	s1 = parsed_args["s1"]
	day0 = parsed_args["day0"]
	day1 = parsed_args["day1"]


	outageScheduleAssessment(path, filesch, s0, s1, day0, day1, M)

end

main()
