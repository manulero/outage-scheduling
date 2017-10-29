
@everywhere function readData(filename)

    data, headers = readcsv(filename; header=true)
    Data = Dict()
    for i in eachindex(headers)
        push!(Data, Symbol(headers[i]) => data[:,i])
    end

    return Data

end


@everywhere function readNetworkData(path)


    NetworkData = Dict()
    push!(NetworkData, :BusData => readData("$(path)/BusData.csv"))
    push!(NetworkData, :LoadData => readData("$(path)/LoadData.csv"))
    push!(NetworkData, :GeneratorData => readData("$(path)/GeneratorData.csv"))
    push!(NetworkData, :GenLookupData => readData("$(path)/GenLookupData.csv"))
    push!(NetworkData, :BranchData => readData("$(path)/BranchData.csv"))
    push!(NetworkData, :WeeklyPeakLoadData => readData("$(path)/WeeklyPeakLoadData.csv"))
    push!(NetworkData, :DailyPeakLoadData => readData("$(path)/DailyPeakLoadData.csv"))
    push!(NetworkData, :HourlyPeakLoadData => readData("$(path)/HourlyPeakLoadData.csv"))

    return NetworkData

end


@everywhere function writeMicroData(path, s, day0, day1, outG, HCap, Ploadfc, Qloadfc, Pload, Qload, uMC, PMC, vMC)


    println("Writing data of micro-scenario $s from days $(day0) to $(day1)")

    outdir = string(path, "/micro-scenarios")
    if !isdir(outdir)
        mkdir(outdir)
    end

    outdirmicro = string(outdir, "/micro$(s)")
    if !isdir(outdirmicro)
        mkdir(outdirmicro)
    end

    for day = day0:day1
        outdirday = string(outdirmicro, "/day$(day)")
        if !isdir(outdirday)
            mkdir(outdirday)
        end

        outfile = string(outdirday, "/outG.csv")
        writecsv(outfile, outG[day])

        outfile = string(outdirday, "/HCap.csv")
        writecsv(outfile, HCap[day])

        outfile = string(outdirday, "/Ploadfc.csv")
        writecsv(outfile, round.(Ploadfc[day], 3))

        outfile = string(outdirday, "/Qloadfc.csv")
        writecsv(outfile, round.(Qloadfc[day], 3))

        outfile = string(outdirday, "/Pload.csv")
        writecsv(outfile, round.(Pload[day], 3))
        
        outfile = string(outdirday, "/Qload.csv")
        writecsv(outfile, round.(Qload[day], 3))

        outfile = string(outdirday, "/uMC.csv")
        I, J, V = findnz(uMC[day])
        if isempty(I)
            I = [1]
            J = [1]
            V = [0]
        end
        writecsv(outfile, (I, J, V))   

        outfile = string(outdirday, "/PMC.csv")
        I, J, V = findnz(PMC[day])
        if isempty(I)
            I = [1]
            J = [1]
            V = [0.0]
        end 
        V = round.(V, 3)
        writecsv(outfile, (I, J, V))

        outfile = string(outdirday, "/vMC.csv")
        I, J, V = findnz(vMC[day])
        if isempty(I)
            I = [1]
            J = [1]
            V = [0]
        end
        writecsv(outfile, (I, J, V))
    end

end



@everywhere function readMicroData(path, s, day)


    MicroData = Dict()

    indir = string(path, "/micro-scenarios/micro$(s)/day$(day)")
    if(!isdir(indir))
        return error("directory $(indir) not found!")
    end

    infile = string(indir, "/outG.csv")
    if(!isfile(infile))
        return error("file $(infile) not found!")
    end
    outG = try
        readcsv(infile, Int)
    catch
        Int[]
    end

    infile = string(indir, "/HCap.csv")
    if(!isfile(infile)) 
        return error("file $(infile) not found!")
    end
    HCap = readcsv(infile, Float64)

    infile = string(indir, "/Ploadfc.csv")
    if(!isfile(infile)) 
        return error("file $(infile) not found!")
    end
    Ploadfc = readcsv(infile, Float64)

    infile = string(indir, "/Qloadfc.csv")
    if(!isfile(infile)) 
        return error("file $(infile) not found!")
    end
    Qloadfc = readcsv(infile, Float64)

    infile = string(indir, "/Pload.csv")
    if(!isfile(infile)) 
        return error("file $(infile) not found!")
    end
    Pload = readcsv(infile, Float64)

    infile = string(indir, "/Qload.csv")
    if(!isfile(infile)) 
        return error("file $(infile) not found!")
    end
    Qload = readcsv(infile, Float64)

    infile = string(indir, "/uMC.csv")
    if(!isfile(infile)) 
        return error("file $(infile) not found!")
    end
    uMCdata = readcsv(infile)
    I = round.(Int64, vec(uMCdata[1,:]))
    J = round.(Int64, vec(uMCdata[2,:]))
    V = ones(Int64, length(I))
    uMCsparse = sparse(I,J,V)

    infile = string(indir, "/PMC.csv")
    if(!isfile(infile)) 
        return error("file $(infile) not found!")
    end    
    PMCdata = readcsv(infile)
    I = round.(Int64, vec(PMCdata[1,:]))
    J = round.(Int64, vec(PMCdata[2,:]))
    V = vec(PMCdata[3,:])
    PMCsparse = sparse(I,J,V)

    infile = string(indir, "/vMC.csv")
    if(!isfile(infile)) 
        return error("file $(infile) not found!")
    end 
    vMCdata = readcsv(infile)
    I = round.(Int64, vec(vMCdata[1,:]))
    J = round.(Int64, vec(vMCdata[2,:]))
    V = ones(Int64, length(I))
    vMCsparse = sparse(I,J,V)

    push!(MicroData, :outG => outG)
    push!(MicroData, :HCap => HCap)
    push!(MicroData, :Ploadfc => Ploadfc)
    push!(MicroData, :Qloadfc => Qloadfc)
    push!(MicroData, :Pload => Pload)
    push!(MicroData, :Qload => Qload)
    push!(MicroData, :uMCsparse => uMCsparse)
    push!(MicroData, :PMCsparse => PMCsparse)
    push!(MicroData, :vMCsparse => vMCsparse)

    return MicroData

end


@everywhere function writeDAoutcome(path, s, day, outB, DAoutcome)

    outBstr = "outB"
    for i=1:length(outB)
        outBstr = string(outBstr, "_$(outB[i])")
    end
    println("Writing $(outBstr) DA data of micro-scenario $s at day $(day)")

    outdirmicro = string(path, "/micro-scenarios/micro$(s)")

    if !isdir(outdirmicro)
        mkdir(outdirmicro)
    end

    outdirday = string(outdirmicro, "/day$(day)")
    if !isdir(outdirday)
        mkdir(outdirday)
    end

    outdiroutB = string(outdirday, "/", outBstr)
    if !isdir(outdiroutB)
        mkdir(outdiroutB)
    end

    outfile = string(outdiroutB, "/DAcost.csv")
    writecsv(outfile, round(DAoutcome[:DAcost], 3))

    outfile = string(outdiroutB, "/uDA.csv")
    I, J, V = findnz(DAoutcome[:uDAsparse])
    if isempty(I)
        I = [1]
        J = [1]
        V = [0]
    end
    writecsv(outfile, (I, J, V)) 

    outfile = string(outdiroutB, "/PDA.csv")
    I, J, V = findnz(DAoutcome[:PDAsparse])
    if isempty(I)
        I = [1]
        J = [1]
        V = [0.0]
    end
    V = round.(V, 3)
    writecsv(outfile, (I, J, V)) 

    outfile = string(outdiroutB, "/DAexec_time.csv")
    writecsv(outfile, round(DAoutcome[:DAexec_time], 3))   

end


@everywhere function readDAoutcome(path, s, day, outB)

    DAoutcome = Dict()

    outBstr = "outB"
    for i=1:length(outB)
        outBstr = string(outBstr, "_$(outB[i])")
    end

    indir = string(path, "/micro-scenarios/micro$(s)/day$(day)/$(outBstr)")

    infile = string(indir, "/DAcost.csv")
    DAcost = try 
        readcsv(infile, Float64)
    catch err
        println(err)
        nothing
    end

    infile = string(indir, "/uDA.csv")
    uDAsparse = try
        uDAdata = readcsv(infile)
        I = round.(Int64, vec(uDAdata[1,:]))
        J = round.(Int64, vec(uDAdata[2,:]))
        V = ones(Int64, length(I))
        sparse(I,J,V)
    catch err
        println(err)
        nothing
    end

    infile = string(indir, "/PDA.csv")
    PDAsparse = try
        PDAdata = readcsv(infile)
        I = round.(Int64, vec(PDAdata[1,:]))
        J = round.(Int64, vec(PDAdata[2,:]))
        V = vec(PDAdata[3,:])
        sparse(I,J,V)
    catch err
        println(err)
        nothing
    end

    push!(DAoutcome, :DAcost => DAcost)
    push!(DAoutcome, :uDAsparse => uDAsparse)
    push!(DAoutcome, :PDAsparse => PDAsparse)


    return DAoutcome

end


@everywhere function writeRToutcome(path, s, day, outB, RToutcome)

    outBstr = "outB"
    for i=1:length(outB)
        outBstr = string(outBstr, "_$(outB[i])")
    end
    println("Writing $(outBstr) data of micro-scenario $s at day $(day)")

    outdirmicro = string(path, "/micro-scenarios/micro$(s)")

    if !isdir(outdirmicro)
        mkdir(outdirmicro)
    end

    outdirday = string(outdirmicro, "/day$(day)")
    if !isdir(outdirday)
        mkdir(outdirday)
    end

    outdiroutB = string(outdirday, "/", outBstr)
    if !isdir(outdiroutB)
        mkdir(outdiroutB)
    end

    outfile = string(outdiroutB, "/RTcost.csv")
    writecsv(outfile, round.(RToutcome[:RTcost], 3))

    outfile = string(outdiroutB, "/PRTdiff0.csv")
    I, J, V = findnz(RToutcome[:PRTdiff0sparse])
    if isempty(I)
        I = [1]
        J = [1]
        V = [0.0]
    end
    V = round.(V, 3)
    writecsv(outfile, (I, J, V)) 

    PRTdiffcsparsedict = RToutcome[:PRTdiffcsparsedict]
    for (key, value) in PRTdiffcsparsedict
        outfile = string(outdiroutB, "/PRTdiffc$(key).csv")
        I, J, V = findnz(value)
        if isempty(I)
            I = [1]
            J = [1]
            V = [0.0]
        end    
        V = round.(V, 3)
        writecsv(outfile, (I, J, V))     
    end


    outfile = string(outdiroutB, "/PRTShed0.csv")
    I, J, V = findnz(RToutcome[:PRTShed0sparse])
    if isempty(I)
        I = [1]
        J = [1]
        V = [0.0]
    end
    V = round.(V, 3)
    writecsv(outfile, (I, J, V)) 


    PRTShedcsparsedict = RToutcome[:PRTShedcsparsedict]
    for (key, value) in PRTShedcsparsedict
        outfile = string(outdiroutB, "/PRTShedc$(key).csv")
        I, J, V = findnz(value)
        if isempty(I)
            I = [1]
            J = [1]
            V = [0.0]
        end    
        V = round.(V, 3)
        writecsv(outfile, (I, J, V))     
    end

    outfile = string(outdiroutB, "/RTexec_time.csv")
    writecsv(outfile, round.(RToutcome[:RTexec_time], 3))
    

end


@everywhere function readRToutcome(path, s, day, outB)

    RToutcome = Dict()

    outBstr = "outB"
    for i=1:length(outB)
        outBstr = string(outBstr, "_$(outB[i])")
    end

    indir = string(path, "/micro-scenarios/micro$(s)/day$(day)/$(outBstr)")

    infile = string(indir, "/RTcost.csv")
    RTcost = try 
        readcsv(infile, Float64)
    catch err
        println(err)
        nothing
    end

    infile = string(indir, "/PRTdiff0.csv")
    PRTdiff0sparse = try
        PRTdiff0data = readcsv(infile)
        I = round.(Int64, vec(PRTdiff0data[1,:]))
        J = round.(Int64, vec(PRTdiff0data[2,:]))
        V = vec(PRTdiff0data[3,:])
        sparse(I,J,V)
    catch err
        println(err)
        nothing
    end

    PRTdiffcsparsedict = Dict()
    list = try 
        filter(x->contains(x,"PRTdiffc"), readdir(indir)) 
    catch err
        println(err)
        []
    end

    for filename in list
        infile = string(indir, "/", filename)
        PRTdiffcdata = readcsv(infile)
        I = round.(Int64, vec(PRTdiffcdata[1,:]))
        J = round.(Int64, vec(PRTdiffcdata[2,:]))
        V = vec(PRTdiffcdata[3,:])
        PRTdiffcsparse = sparse(I,J,V)
        c = parse(Int, filename[9:end-4]) #the 9 is the (length of "PRTdiffc") + 1 ; the 4 is the length of ".csv"
        push!(PRTdiffcsparsedict, c => PRTdiffcsparse)
    end

    if isempty(PRTdiffcsparsedict) PRTdiffcsparsedict = nothing end


    infile = string(indir, "/PRTShed0.csv")
    PRTShed0sparse = try
        PRTShed0data = readcsv(infile)
        I = round.(Int64, vec(PRTShed0data[1,:]))
        J = round.(Int64, vec(PRTShed0data[2,:]))
        V = vec(PRTShed0data[3,:])
        sparse(I,J,V)
    catch err
        println(err)
        nothing
    end  

    PRTShedcsparsedict = Dict()
    list = try 
        filter(x->contains(x,"PRTShedc"), readdir(indir)) 
    catch err
        println(err)
        []
    end

    for filename in list
        infile = string(indir, "/", filename)
        PRTShedcdata = readcsv(infile)
        I = round.(Int64, vec(PRTShedcdata[1,:]))
        J = round.(Int64, vec(PRTShedcdata[2,:]))
        V = vec(PRTShedcdata[3,:])
        PRTShedcsparse = sparse(I,J,V)
        c = parse(Int, filename[9:end-4]) #the 9 is the (length of "PRTdiffc") + 1 ; the 4 is the length of ".csv"
        push!(PRTShedcsparsedict, c => PRTShedcsparse)
    end
    if isempty(PRTShedcsparsedict) PRTShedcsparsedict = nothing end

    push!(RToutcome, :RTcost => RTcost)
    push!(RToutcome, :PRTdiff0sparse => PRTdiff0sparse)
    push!(RToutcome, :PRTdiffcsparsedict => PRTdiffcsparsedict)
    push!(RToutcome, :PRTShed0sparse => PRTShed0sparse)
    push!(RToutcome, :PRTShedcsparsedict => PRTShedcsparsedict)

    return RToutcome

end

@everywhere function writeAssessmentData(path, s, day, outB, assessment, assessment_time, kwrd)

    outBstr = "outB"
    for i=1:length(outB)
        outBstr = string(outBstr, "_$(outB[i])")
    end
    println("Writing $(outBstr) assessment data of micro-scenario $s at day $(day)")  

    outdirmicro = string(path, "/micro-scenarios/micro$(s)")

    if !isdir(outdirmicro)
        mkdir(outdirmicro)
    end

    outdirday = string(outdirmicro, "/day$(day)")
    if !isdir(outdirday)
        mkdir(outdirday)
    end

    outdiroutB = string(outdirday, "/", outBstr)
    if !isdir(outdiroutB)
        mkdir(outdiroutB)
    end

    outfile = string(outdiroutB, "/$(kwrd)_assessment.csv")
    writecsv(outfile, round.(assessment, 3))

    outfile = string(outdiroutB, "/$(kwrd)_assessment_time.csv")
    writecsv(outfile, round.(assessment_time, 3))
    

end


@everywhere function readAssessmentData(path, s, day, outB, kwrd)

    outBstr = "outB"
    for i=1:length(outB)
        outBstr = string(outBstr, "_$(outB[i])")
    end

    indir = string(path, "/micro-scenarios/micro$(s)/day$(day)/$(outBstr)")

    infile = string(indir, "/$(kwrd)_assessment.csv")
    assessment = readcsv(infile, Float64)

    return assessment

end



function readOutageCostData(path, s, day, outB)

    outBstr = "outB"
    for i=1:length(outB)
        outBstr = string(outBstr, "_$(outB[i])")
    end

    indir = string(path, "/micro-scenarios/micro$(s)/day$(day)/$(outBstr)")

    infileDA = string(indir, "/DAcost.csv")
    infileRT = string(indir, "/RTcost.csv")
    DAcost = readcsv(infileDA, Float64)[1]
    RTcost = readcsv(infileRT, Float64)

    return DAcost, RTcost

end 

function readOutageAssessmentData(path, s, day, outB)

    outBstr = "outB"
    for i=1:length(outB)
        outBstr = string(outBstr, "_$(outB[i])")
    end

    indir = string(path, "/micro-scenarios/micro$(s)/day$(day)/$(outBstr)")

    DAcost, RTcost = readOutageCostData(path, s, day, outB)

    infile_assessment = string(indir, "/assessment_cost.csv")
    assessment_cost = readcsv(infile_assessment, Float64)

    return DAcost, RTcost, assessment_cost

end



#=
@everywhere function readBranchDataOld(inputfile, N)


    println("Reading branch data from $(inputfile)")

    local 
        beta = [],
        X = [],
        fmax = [],
        STfmax = [],
        BranchOutRate = [],
        BranchOutDuration = [],
        L
    open(inputfile) do filehandle

        key = "BRANCH DATA"
        len = length(key)
        while !eof(filehandle)
            line = readline(filehandle)
            if(length(line) >= len && line[1:len] == key)
                break
            end
        end
        while !eof(filehandle)
            line = readline(filehandle)
            if(length(line) >= 3 && line[1:3] == "END")
                break
            end
            splitted = split(line," ")
            temp = zeros(N)
            fr_bus = parse(Int, splitted[1])
            to_bus = parse(Int, splitted[2])
            temp[fr_bus] = 1
            temp[to_bus] = -1
            append!(beta, temp)
            push!(X, float(splitted[3]))
            push!(fmax, float(splitted[4]))
            push!(STfmax, float(splitted[5]))
            push!(BranchOutRate, float(splitted[6]))
            push!(BranchOutDuration, parse(Int, splitted[7]))
        end
        L = length(X)
        beta = reshape(beta, N, L)

    end
    
    return beta, X, fmax, STfmax, BranchOutRate, BranchOutDuration

end

#######Micro-scenario data

@everywhere function readGenAvailability(inputfile, G)

    println("Reading generator availability from $(inputfile)")

    local GenSchOut = Array(Any, G)
    local GenForceOut = Array(Any, G)
    for g=1:G
        GenSchOut[g] = Int[]
        GenForceOut[g] = Tuple{Int, Int}[]
    end
    
    open(inputfile) do filehandle

        key = "GENERATOR SCHEDULED OUTAGES"
        len = length(key)
        while !eof(filehandle)
            line = readline(filehandle)
            if(length(line) >= len && line[1:len] == key)
                break
            end
        end
        while !eof(filehandle)
            line = readline(filehandle)
            if(length(line) >= 3 && line[1:3] == "END")
                break
            end
            splitted = split(line," ")
            g = parse(Int64, splitted[1])
            append!(GenSchOut[g], [parse(Int64,s) for s = splitted[2:end]])
        end
    end

    open(inputfile) do filehandle
        key = "GENERATOR FORCED OUTAGES"
        len = length(key)
        while !eof(filehandle)
            line = readline(filehandle)
            if(length(line) >= len && line[1:len] == key)
                break
            end
        end
        while !eof(filehandle)
            line = readline(filehandle)
            if(length(line) >= 3 && line[1:3] == "END")
                break
            end
            splitted = split(line," ")
            g = parse(Int64, splitted[1])
            xx = [split(strip(chomp(s),['(',')']),",") for s = splitted[2:end]]
            append!(GenForceOut[g], Tuple{Int,Int}[tuple(parse(Int64, yy[1]), parse(Int64, yy[2])) for yy = xx])
        end
    end
    
    return GenSchOut, GenForceOut
    
end

@everywhere function readBranchAvailability(inputfile, L)

    println("Reading branch availability from $(inputfile)")

    #local BranchSchOut = Array(Any, L)
    local BranchForceOut = Array(Any, L)
    for l=1:L
        BranchForceOut[l] = Tuple{Int, Int}[]
    end
    open(inputfile) do filehandle

        key = "BRANCH FORCED OUTAGES"
        len = length(key)
        while !eof(filehandle)
            line = readline(filehandle)
            if(length(line) >= len && line[1:len] == key)
                break
            end
        end
        while !eof(filehandle)
            line = readline(filehandle)
            if(length(line) >= 3 && line[1:3] == "END")
                break
            end
            splitted = split(line," ")
            l = parse(Int64, splitted[1])
            xx = [split(strip(chomp(s),['(',')']),",") for s = splitted[2:end]]
            append!(BranchForceOut[l], Tuple{Int,Int}[tuple(parse(Int64, yy[1]), parse(Int64, yy[2])) for yy = xx])

        end
    end
    return BranchForceOut
end




@everywhere function readDemandRealization(inputfile)

    println("Reading demand realization from $(inputfile)")

    local Pload = []     #pu

    open(inputfile) do filehandle

        key = "DEMAND REALIZATION"
        len = length(key)
        while !eof(filehandle)
            line = readline(filehandle)
            if(length(line) >= len && line[1:len] == key)
                break
            end
        end
        i = 1
        while !eof(filehandle)
            line = readline(filehandle)
            if(length(line) >= 3 && line[1:3] == "END")
                break
            end
            splitted = split(line," ")
            append!(Pload, float(splitted))     #append means row vector
            i = i + 1
        end
        D = i - 1
        H = round(Int, length(Pload) / D)
        Pload = reshape(Pload, H, D)
        Pload = transpose(Pload)

    end

    return Pload

end



@everywhere function readDemandForecast(inputfile)

    println("Reading demand forecast from $(inputfile)")

    local Ploadfc = []     #pu

    open(inputfile) do filehandle

        key = "DEMAND FORECAST"
        len = length(key)
        while !eof(filehandle)
            line = readline(filehandle)
            if(length(line) >= len && line[1:len] == key)
                break
            end
        end
        i = 1
        while !eof(filehandle)
            line = readline(filehandle)
            if(length(line) >= 3 && line[1:3] == "END")
                break
            end
            splitted = split(line," ")
            append!(Ploadfc, float(splitted))
            i = i + 1
        end
        D = i - 1
        H = round(Int, length(Ploadfc) / D)
        Ploadfc = reshape(Ploadfc, H, D)
        Ploadfc = transpose(Ploadfc)

    end

    return Ploadfc

end



@everywhere function ReadHydroCapacity(inputfile)

    println("Reading hydro capacity from $(inputfile)")

    local HHourCap = []     #pu

    open(inputfile) do filehandle

        key = "HYDRO CAPACITY"
        len = length(key)
        while !eof(filehandle)
            line = readline(filehandle)
            if(length(line) >= len && line[1:len] == key)
                break
            end
        end
        while !eof(filehandle)
            line = readline(filehandle)
            if(length(line) >= 3 && line[1:3] == "END")
                break
            end
            splitted = split(line," ")
            append!(HHourCap, float(splitted))
        end

    end

    return HHourCap

end


@everywhere function readMarketOutcome(inputfile, G)

    println("Reading market outcome from $(inputfile)")

    local MarketClearing = Array(Any, G)
    for g=1:G
        MarketClearing[g] = Tuple{Int, Int, Float64}[]
    end

    open(inputfile) do filehandle

        key = "MARKET CLEARING OUTCOME"
        len = length(key)
        while !eof(filehandle)
            line = readline(filehandle)
            if(length(line) >= len && line[1:len] == key)
                break
            end
        end
        while !eof(filehandle)
            line = readline(filehandle)
            if(length(line) >= 3 && line[1:3] == "END")
                break
            end
            splitted = split(line," ")
            g = parse(Int64, splitted[1])
            xx = [split(strip(chomp(s),['(',')']),",") for s = splitted[2:end]]
            append!(MarketClearing[g], Tuple{Int,Int,Float64}[tuple(parse(Int64, yy[1]), parse(Int64, yy[2]), parse(Float64, yy[3])) for yy = xx])
        end

    end
    return MarketClearing
end

=#

#=

function readMaintPolicyData(inputfile)

    local
        LineOutages = []
        LineOutDuration = []

    open(inputfile) do filehandle

        key = "MAINTENANCE POLICY DATA"
        len = length(key)
        while !eof(filehandle)
            line = readline(filehandle)
            if(length(line) >= len && line[1:len] == key)
                break
            end
        end
        while !eof(filehandle)
            line = readline(filehandle)
            if(length(line) >= 3 && line[1:3] == "END")
                break
            end
            splitted = split(line," ")
            push!(LineOutages, parse(Int64, splitted[1]))
            push!(LineOutDuration, parse(Int64, chomp(splitted[2])))
        end

    end
    return LineOutages, LineOutDuration
end

function readOutageCostData2(inputfile, nb_scenarios)

    println("Reading outage cost data from $(inputfile) with $(nb_scenarios) scenarios")


    hpd = 24
    local
        DADailyCostj, RTHourlyCostj


    open(inputfile) do filehandle

        key = "Days"
        len = length(key)
        DAkey = "DADailyCost"
        DAlen = length(DAkey)
        RTkey = "RTHourlyCost"
        RTlen = length(RTkey)
        DAdone = 0
        RTdone = 0
        while !eof(filehandle)
            if DAdone == nb_scenarios && RTdone == nb_scenarios
                break
            end
            line = readline(filehandle)
            splitted = split(line, " ")
            if(length(line) >= len && line[1:len] == key)
                days = parse(Int, splitted[4]) - parse(Int, splitted[2]) + 1
                DADailyCostj = zeros(days)
                RTHourlyCostj = zeros(days * hpd)
            end
            if length(splitted) < 5 continue end
            if(length(splitted[5]) >= DAlen && splitted[5][1:DAlen] == DAkey && DAdone < nb_scenarios)
                splitted2 = split(splitted[5], [',', '[', ']'])
                DADailyCostj += float(splitted2[2:end-1])
                DAdone += 1
            end
            if(length(splitted[5]) >= RTlen && splitted[5][1:RTlen] == RTkey && RTdone < nb_scenarios)
                splitted2 = split(splitted[5], [',', '[', ']'])
                RTHourlyCostj += float(splitted2[2:end-1])
                RTdone += 1
            end
        end
    end
    days = round(Int, length(RTHourlyCostj) / hpd)
    RTDailyCostj = Array{Float64}(days)
    for d=1:days
        RTDailyCostj[d] = sum(RTHourlyCostj[(d-1)*hpd+1:d*hpd])
    end


    return (DADailyCostj + RTDailyCostj) / nb_scenarios

end

# scenario_offset tells how many lines do I have to read before entering the information
# day_offset does a similar thing
function readOutageCostData(inputfile, day0, day1, s0, s1)

    println("Reading outage cost data from $(inputfile)")

    splitted = split(inputfile, "_")

    f_days = split(splitted[2][5:end], "-")
    f_day0 = parse(Int, f_days[1])
    f_day1 = parse(Int, f_days[2])

    f_scenarios = split(splitted[3][10:end], "-")
    f_s0 = parse(Int, f_scenarios[1])
    f_s1 = parse(Int, f_scenarios[2])

    #Base.run(`mac2unix $(inputfile)`)

    local DADailyCost = [], RTHourlyCost = []

    open(inputfile) do filehandle

        key = "DA daily cost data"
        len = length(key)
        while !eof(filehandle)
            line = readline(filehandle)
            if(length(line) >= len && line[1:len] == key)
                break
            end
        end
        s = f_s0
        while !eof(filehandle)
            if (s > s1)
                break
            end
            line = readline(filehandle)
            if (length(line) >= 3 && line[1:3] == "END")
                break
            end
            if (s >= s0)
                splitted = split(line, " ")
                push!(DADailyCost, float(splitted))
                s += 1
            end
        end

        key = "RT hourly cost data"
        len = length(key)
        while !eof(filehandle)
            line = readline(filehandle)
            if(length(line) >= len && line[1:len] == key)
                break
            end
        end
        while !eof(filehandle)
            line = readline(filehandle)
            if(length(line) >= 3 && line[1:3] == "END")
                break
            end
            splitted = split(line, " ")
            push!(RTHourlyCost, float(splitted))
        end
    end

    return DADailyCost, RTHourlyCost

end


#helper function to read all network data
@everywhere function readNetworkData(fn)

    #read non-stochastic data
    N = readNumberOfBuses(fn)
    Dn, D, PCTBusLoad, VoLL = readLoadData(fn, N)
    Gn, Gtype, G = readGeneratorData(fn, N)
    Plim, RD, RU, DT, UT, su_cost, cost = readGenLookupTable(fn)
    beta, X, fmax, STfmax, L, BranchOutRate, BranchOutDuration = readBranchData(fn, N)
    PCTWeeklyPeakLoad, PCTDailyPeakLoad, PCTHourlyPeakLoad = readPeakLoadData(fn)


    return N, Dn, D, Gn, Gtype, G, Plim, RD, RU, DT, UT, su_cost, cost, beta, X, fmax, STfmax, L, BranchOutRate, BranchOutDuration, PCTBusLoad, PCTWeeklyPeakLoad, PCTDailyPeakLoad, PCTHourlyPeakLoad, VoLL
end

#helper function to read all micro-scenario data
@everywhere function readMicroScenarioData(fn, G, L)

    GenSchOut, GenForceOut = readGenAvailability(fn, G)
    BranchForceOut = readBranchAvailability(fn, L)
    MarketClearing = readMarketOutcome(fn, G)
    Ploadfc = readDemandForecast(fn)
    Pload = readDemandRealization(fn)
    HHydroCap = ReadHydroCapacity(fn)

    return GenSchOut, GenForceOut, Ploadfc, MarketClearing, BranchForceOut, HHydroCap, Pload
end


=#

#=


@everywhere function readMicroData(path, s, day, DataKeys)

    MicroData = Dict()

    indir = string(path, "/micro-scenarios/micro$(s)/day$(day)")
    if(!isdir(indir))
        return error("directory $(indir) not found!")
    end

    for datakey in DataKeys

        infile = string(indir, "/$(datakey).csv")
        if(!isfile(infile))
            return error("file $(infile) not found!")
        end
        data = try
            readcsv(infile)
        catch
            []
        end
        push!(MicroData, datakey => data)
    end
end

=#