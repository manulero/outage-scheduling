#!~/Downloads/julia-0.4.5/usr/bin/julia

using ArgParse
using OffsetArrays

include("./dataInputOutput.jl")
include("./microGenerativeModel.jl")


function generateMicroscenarios(day0, day1, s0, s1, M, path)

    addprocs(M)
    M = nworkers()

    #these includes come after addprocs() so files are loaded on all processes
    include("./dataInputOutput.jl")
    include("./microGenerativeModel.jl")

    ref = OffsetArray(Any, s0:s1)

    for s=s0:s1
        ref[s] = @spawn genMicros(s, day0, day1, path)
    end

    for s=s0:s1
        fetch(ref[s])
    end

    for i in workers()
        rmprocs(i)
    end


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
            default = 1
        "s0"
            help = "The first micro-scenario"
            arg_type = Int
            default = 1
        "s1"
            help = "The last micro-scenario"
            arg_type = Int
            default = 4
        "day0"
            help = "The first day"
            arg_type = Int
            default = 1
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
    s0 = parsed_args["s0"]
    s1 = parsed_args["s1"]
    day0 = parsed_args["day0"]
    day1 = parsed_args["day1"]


    generateMicroscenarios(day0, day1, s0, s1, M, path)

end

main()
