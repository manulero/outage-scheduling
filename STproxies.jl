#!~/Codes/julia-0.6.0/bin/julia


@everywhere using JuMP
@everywhere using CPLEX
@everywhere using Ipopt
@everywhere using OffsetArrays

@everywhere include("./dataInputOutput.jl")

@everywhere function safesplit(data)

    return try split(data, " ") catch; data end

end

@everywhere function fillWithZeros(A, nrows, ncols)

    return [full(A) zeros(eltype(A), size(A,1), ncols - size(A,2)); zeros(eltype(A), nrows - size(A,1), ncols)]

end



@everywhere function prepareOPFNetworkData(NetworkData)


    OPFNetworkData = Dict()

    #first split NetworkData into components
    BusData = NetworkData[:BusData]
    LoadData = NetworkData[:LoadData]
    GeneratorData = NetworkData[:GeneratorData]
    GenLookupData = NetworkData[:GenLookupData]
    BranchData = NetworkData[:BranchData]
    WeeklyPeakLoadData = NetworkData[:WeeklyPeakLoadData]
    DailyPeakLoadData = NetworkData[:DailyPeakLoadData]
    HourlyPeakLoadData = NetworkData[:HourlyPeakLoadData]

    #now fill OPFNetworkData from NetworkData components
    Btype = BusData[:Type]
    push!(OPFNetworkData, :Btype => Btype)
    push!(OPFNetworkData, :Vsp => BusData[:Vsp])
    push!(OPFNetworkData, :Bgs => BusData[:gs])
    push!(OPFNetworkData, :Bbs => BusData[:bs])
    Gbusid = round.(Int, GeneratorData[:Busid])
    push!(OPFNetworkData, :Gbusid => Gbusid)
    push!(OPFNetworkData, :Gtype => GeneratorData[:Type])
    push!(OPFNetworkData, :RampDn => Dict(zip(GenLookupData[:Type], GenLookupData[:RampDn])))
    push!(OPFNetworkData, :RampUp => Dict(zip(GenLookupData[:Type], GenLookupData[:RampDn])))
    push!(OPFNetworkData, :MinDnTime => Dict(zip(GenLookupData[:Type], GenLookupData[:MinDnTime])))
    push!(OPFNetworkData, :MinUpTime => Dict(zip(GenLookupData[:Type], GenLookupData[:MinUpTime])))
    push!(OPFNetworkData, :SUCost => Dict(zip(GenLookupData[:Type], GenLookupData[:fuelPrice] .* GenLookupData[:SUCost])))
    push!(OPFNetworkData, :fuelType => Dict(zip(GenLookupData[:Type], GenLookupData[:fuelType])))
    push!(OPFNetworkData, :Pout => Dict(zip(GenLookupData[:Type], [float(safesplit(GenLookupData[:Pout][i])) for i=1:length(GenLookupData[:Pout])])))
    push!(OPFNetworkData, :Pmin => Dict(zip(GenLookupData[:Type], [float(safesplit(GenLookupData[:Pout][i]))[1] for i=1:length(GenLookupData[:Pout])])))
    Pmax = Dict(zip(GenLookupData[:Type], [sum(float(safesplit(GenLookupData[:Pout][i]))) for i=1:length(GenLookupData[:Pout])]))
    push!(OPFNetworkData, :Pmax => Pmax)
    push!(OPFNetworkData, :RampSD => Pmax)  #assume startup and shutdown ramp rates are equal to the Pmax
    push!(OPFNetworkData, :RampSU => Pmax)
    push!(OPFNetworkData, :Qmin => Dict(zip(GenLookupData[:Type], GenLookupData[:Qmin])))
    push!(OPFNetworkData, :Qmax => Dict(zip(GenLookupData[:Type], GenLookupData[:Qmax])))
    push!(OPFNetworkData, :Cost => Dict(zip(GenLookupData[:Type], GenLookupData[:fuelPrice] .* [float(safesplit(GenLookupData[:Cost][i])) for i=1:length(GenLookupData[:Cost])])))
    Dbusid = round.(Int, LoadData[:Busid])
    push!(OPFNetworkData, :Dbusid => Dbusid)
    push!(OPFNetworkData, :VoLL => LoadData[:VoLL])
    FrBusid = round.(Int, BranchData[:FrBusid])
    ToBusid = round.(Int, BranchData[:ToBusid])
    push!(OPFNetworkData, :FrBusid => FrBusid)
    push!(OPFNetworkData, :ToBusid => ToBusid)
    push!(OPFNetworkData, :r => BranchData[:r])
    push!(OPFNetworkData, :x => BranchData[:x])
    push!(OPFNetworkData, :b => BranchData[:b])
    push!(OPFNetworkData, :tr => BranchData[:tr])
    push!(OPFNetworkData, :fmax => BranchData[:fmax])
    push!(OPFNetworkData, :STfmax => BranchData[:STfmax])
    push!(OPFNetworkData, :BranchOutRate => BranchData[:BranchOutRate])
    N = length(Btype)
    G = length(Gbusid)
    D = length(Dbusid)
    L = length(FrBusid)
    G_n = [Set() for i = 1:N]
    D_n = [Set() for i = 1:N]
    Fr_n = [Set() for i = 1:N]
    To_n = [Set() for i = 1:N]
    for g=1:G push!(G_n[Gbusid[g]], g) end
    for d=1:D push!(D_n[Dbusid[d]], d) end
    for l=1:L push!(Fr_n[FrBusid[l]], l) end
    for l=1:L push!(To_n[ToBusid[l]], l) end
    push!(OPFNetworkData, :N => N)
    push!(OPFNetworkData, :G => G)
    push!(OPFNetworkData, :D => D)
    push!(OPFNetworkData, :L => L)
    push!(OPFNetworkData, :G_n => G_n)
    push!(OPFNetworkData, :D_n => D_n)
    push!(OPFNetworkData, :Fr_n => Fr_n)
    push!(OPFNetworkData, :To_n => To_n)
    push!(OPFNetworkData, :hpd => length(HourlyPeakLoadData[:PCTLoadWinterWeekday]))
    push!(OPFNetworkData, :dpw => length(DailyPeakLoadData[:PCTLoad]))
    push!(OPFNetworkData, :wpy => length(WeeklyPeakLoadData[:PCTLoad]))


    return OPFNetworkData

end

@everywhere function DANminus1proxy(NetworkData, MicroData, outB)


    OPFNetworkData = prepareOPFNetworkData(NetworkData)

    Btype = OPFNetworkData[:Btype]
    Gbusid = OPFNetworkData[:Gbusid]
    Gtype = OPFNetworkData[:Gtype]
    RampDn = OPFNetworkData[:RampDn]
    RampUp = OPFNetworkData[:RampUp]
    MinDnTime = OPFNetworkData[:MinDnTime]
    MinUpTime = OPFNetworkData[:MinUpTime]
    SUCost = OPFNetworkData[:SUCost]
    fuelType = OPFNetworkData[:fuelType]
    Pmin = OPFNetworkData[:Pmin]
    Pmax = OPFNetworkData[:Pmax]
    Cost = OPFNetworkData[:Cost]
    RampSD = OPFNetworkData[:RampSD]
    RampSU = OPFNetworkData[:RampSU]
    Dbusid = OPFNetworkData[:Dbusid]
    FrBusid = OPFNetworkData[:FrBusid]
    ToBusid = OPFNetworkData[:ToBusid]
    x = OPFNetworkData[:x]
    fmax = OPFNetworkData[:fmax]
    STfmax = OPFNetworkData[:STfmax]
    N = OPFNetworkData[:N]
    G = OPFNetworkData[:G]
    D = OPFNetworkData[:D]
    L = OPFNetworkData[:L]
    G_n = OPFNetworkData[:G_n]
    D_n = OPFNetworkData[:D_n]
    Fr_n = OPFNetworkData[:Fr_n]
    To_n = OPFNetworkData[:To_n]
    hpd = OPFNetworkData[:hpd]

    outG = MicroData[:outG]
    HCap = MicroData[:HCap]
    Ploadfc = MicroData[:Ploadfc]
    uMC = fillWithZeros(MicroData[:uMCsparse], G, hpd)
    PMC = fillWithZeros(MicroData[:PMCsparse], G, hpd)
    vMC = fillWithZeros(MicroData[:vMCsparse], G, hpd)
    aGen = ones(Int, G)
    aGen[outG] = 0
    for g=1:G, h=1:hpd
        PMC[g,h] = (uMC[g,h] == 1)? min(max(PMC[g,h], Pmin[Gtype[g]]), Pmax[Gtype[g]]) : PMC[g,h] #force PMC within limits 
    end

    #initialize array of branch availability from outage schedule
    aBranch = ones(Int, L)
    aBranch[outB] = 0

    #N - 1 contingencies
    C = L
    a = ones(Float64, L, C) - eye(L, C)

    #initialize gen status and time of 'must keep status'
    uInit = zeros(Int, G)
    MinKeepTime = zeros(Int, G)

    Dt = 60
    dt = 15

    DAoutcome = Dict()

    tic()

    m = Model(solver = CplexSolver(CPX_PARAM_MIPDISPLAY=0, CPX_PARAM_THREADS=16, CPX_PARAM_EPGAP=0.05))

    @variable(m, uDA[1:G,1:hpd], Bin)
    @variable(m, 0 <= vDA[1:G,1:hpd] <= 1)
    @variable(m, 0 <= wDA[1:G,1:hpd] <= 1)
    @variable(m, vSO[1:G,1:hpd] >= 0)
    @variable(m, Pup0[1:G,1:hpd] >= 0)
    @variable(m, Pdn0[1:G,1:hpd] >= 0)
    @variable(m, Pupc[1:G,1:C,1:hpd] >= 0)
    @variable(m, Pdnc[1:G,1:C,1:hpd] >= 0)
    @variable(m, f0[1:L,1:hpd])
    @variable(m, th0[1:N,1:hpd])
    @variable(m, STfc[1:L,1:C,1:hpd])
    @variable(m, STthc[1:N,1:C,1:hpd])
    @variable(m, fc[1:L,1:C,1:hpd])
    @variable(m, thc[1:N,1:C,1:hpd])

    @objective(m, Min, sum(sum(SUCost[Gtype[g]] * vSO[g,h] + Cost[Gtype[g]][end] * Pup0[g,h] * 100 for g=1:G) for h = 1:hpd) )


    @constraint(m, Reference0[h=1:hpd], th0[1,h] == 0 )
    @constraint(m, ReferenceSTc[c=1:C, h=1:hpd], STthc[1,c,h] == 0 )
    @constraint(m, Referencec[c=1:C, h=1:hpd], thc[1,c,h] == 0 )

    #Generator availability and cycling
    @constraint(m, keep[g=1:G, h = 1:MinKeepTime[g]], uDA[g,h] == uInit[g])
    @constraint(m, updown1[g=1:G], vDA[g,1] - wDA[g,1] == uDA[g,1] - uInit[g])
    @constraint(m, updown[g=1:G, h = 2:hpd], vDA[g,h] - wDA[g,h] == uDA[g,h] - uDA[g,h-1] )
    @constraint(m, MinUPTime[g=1:G, h = MinUpTime[Gtype[g]]:hpd], sum(vDA[g,s] for s=h-MinUpTime[Gtype[g]]+1:h) <= uDA[g,h] )
    @constraint(m, MinDNTime[g=1:G, h = MinDnTime[Gtype[g]]:hpd], sum(wDA[g,s] for s=h-MinDnTime[Gtype[g]]+1:h) <= 1 - uDA[g,h] )
    @constraint(m, TSOpay[g=1:G, h = 1:hpd], vSO[g,h] >= vDA[g,h] - vMC[g,h] )


    #Pre-contingency state
    @constraint(m, Gen0Min[g=1:G,h=1:hpd], Pup0[g,h] >= aGen[g] * uDA[g,h] * (Pmin[Gtype[g]] - PMC[g,h]))
    @constraint(m, noHydroGen0Max[g=1:G,h=1:hpd; fuelType[Gtype[g]] != "Hydro"], Pup0[g,h] <= aGen[g] * uDA[g,h] * (Pmax[Gtype[g]] - PMC[g,h]))
    @constraint(m, HydroGen0Max[g=1:G,h=1:hpd; fuelType[Gtype[g]] == "Hydro"], Pup0[g,h] <= aGen[g] * uDA[g,h] * (HCap[h] - PMC[g,h]))
    @constraint(m, Gen0Minbis[g=1:G,h=1:hpd], Pdn0[g,h] <= aGen[g] * uDA[g,h] * uMC[g,h] * (PMC[g,h] - Pmin[Gtype[g]]))


    @constraint(m, Power0[n=1:N, h=1:hpd], 
        sum(uDA[g,h] * PMC[g,h] + Pup0[g,h] - Pdn0[g,h] for g in G_n[n]) - 
        (sum(f0[l,h] for l in Fr_n[n]) - sum(f0[l,h] for l in To_n[n])) == 
        sum(Ploadfc[d,h] for d in D_n[n])
        )
    @constraint(m, RampDownInterPeriod[g=1:G, h=2:hpd], 
        (uDA[g,h-1] * PMC[g,h-1] + Pup0[g,h-1] - Pdn0[g,h-1]) - 
        (uDA[g,h] * PMC[g,h] + Pup0[g,h] - Pdn0[g,h]) <= 
        Dt * RampDn[Gtype[g]] * uDA[g,h] + RampSD[Gtype[g]] * wDA[g,h]
        )
    @constraint(m, RampUpInterPeriod[g=1:G, h=2:hpd], 
        (uDA[g,h] * PMC[g,h] + Pup0[g,h] - Pdn0[g,h]) - 
        (uDA[g,h-1] * PMC[g,h-1] + Pup0[g,h-1] - Pdn0[g,h-1]) <= 
        Dt * RampUp[Gtype[g]] * uDA[g,h-1] + RampSU[Gtype[g]] * vDA[g,h]
        )

    @constraint(m, Flow0[l=1:L, h=1:hpd], f0[l,h] - aBranch[l] * (1/x[l]) * (th0[FrBusid[l],h] - th0[ToBusid[l],h]) == 0 )
    @constraint(m, Flow0PosLim[l=1:L, h=1:hpd], f0[l,h] <= fmax[l] )
    @constraint(m, Flow0NegLim[l=1:L, h=1:hpd], -f0[l,h] <= fmax[l] )


    #Short-term post-contingency state
    @constraint(m, PowerSTc[n=1:N,c=1:C,h=1:hpd], 
        sum(uDA[g,h] * PMC[g,h] + Pup0[g,h] - Pdn0[g,h] for g in G_n[n]) - 
        (sum(STfc[l,c,h] for l in Fr_n[n]) - sum(STfc[l,c,h] for l in To_n[n])) == 
        sum(Ploadfc[d,h] for d in D_n[n])
        )
    @constraint(m, FlowSTc[l=1:L,c=1:C,h=1:hpd], STfc[l,c,h] - aBranch[l] * a[l,c] * (1/x[l]) * (STthc[FrBusid[l],c,h] - STthc[ToBusid[l],c,h]) == 0 )

    @constraint(m, FlowSTcPosLim[l=1:L,c=1:C,h=1:hpd], STfc[l,c,h] <= STfmax[l] )
    @constraint(m, FlowSTcNegLim[l=1:L,c=1:C,h=1:hpd], -STfc[l,c,h] <= STfmax[l] )


    #Post-contingency state (after the successful application of corrective control)
    @constraint(m, Powerc[n=1:N,c=1:C,h=1:hpd], 
        sum(uDA[g,h] * PMC[g,h] + Pup0[g,h] - Pdn0[g,h] + Pupc[g,c,h] - Pdnc[g,c,h] for g in G_n[n]) - 
        (sum(fc[l,c,h] for l in Fr_n[n]) - sum(fc[l,c,h] for l in To_n[n])) == 
        sum(Ploadfc[d,h] for d in D_n[n])
        )
    @constraint(m, Flowc[l=1:L,c=1:C,h=1:hpd], fc[l,c,h] - aBranch[l] * a[l,c] * (1/x[l]) * (thc[FrBusid[l],c,h] - thc[ToBusid[l],c,h]) == 0 )

    @constraint(m, FlowcPosLim[l=1:L,c=1:C,h=1:hpd], fc[l,c,h] <= fmax[l] )
    @constraint(m, FlowcNegLim[l=1:L,c=1:C,h=1:hpd], -fc[l,c,h] <= fmax[l] )
     
    @constraint(m, GencMin[g=1:G,c=1:C,h=1:hpd], Pup0[g,h] - Pdn0[g,h] + Pupc[g,c,h] >= aGen[g] * uDA[g,h] * (Pmin[Gtype[g]] - PMC[g,h]))
    @constraint(m, noHydroGencMax[g=1:G,c=1:C,h=1:hpd; fuelType[Gtype[g]] != "Hydro"], Pup0[g,h] - Pdn0[g,h] + Pupc[g,c,h] <= aGen[g] * uDA[g,h] * (Pmax[Gtype[g]] - PMC[g,h]))
    @constraint(m, HydroGencMax[g=1:G,c=1:C,h=1:hpd; fuelType[Gtype[g]] == "Hydro"], Pup0[g,h] - Pdn0[g,h] + Pupc[g,c,h] <= aGen[g] * uDA[g,h] * (HCap[h] - PMC[g,h]))
    @constraint(m, GencMinbis[g=1:G,c=1:C,h=1:hpd], Pdn0[g,h] - Pup0[g,h] + Pdnc[g,c,h] <= aGen[g] * uDA[g,h] * uMC[g,h] * (PMC[g,h] - Pmin[Gtype[g]]))

    @constraint(m, RampDownc[g=1:G,c=1:C,h=1:hpd], Pdnc[g,c,h] <= dt * RampDn[Gtype[g]] )
    @constraint(m, RampUpc[g=1:G,c=1:C,h=1:hpd], Pupc[g,c,h] <= dt * RampUp[Gtype[g]] )
    
    #TT = STDOUT # save original STDOUT stream
    #redirect_stdout()
    status = solve(m; suppress_warnings=false)
    #redirect_stdout(TT) # restore STDOUT

    if status == :Optimal
        println("solved")
        push!(DAoutcome, :DAcost => getobjectivevalue(m))
        uDAsparse = sparse(round.(Int, getvalue(uDA)))
        push!(DAoutcome, :uDAsparse => uDAsparse)
        push!(DAoutcome, :PDAsparse => sparse(uDAsparse .* PMC + getvalue(Pup0 - Pdn0)))
    else #no solution, just do nothing...
        println("unsolved")
        push!(DAoutcome, :DAcost => 0.0)
        push!(DAoutcome, :uDAsparse => MicroData[:uMCsparse])
        push!(DAoutcome, :PDAsparse => MicroData[:PMCsparse])
    end

    DAexec_time = toq()
    push!(DAoutcome, :DAexec_time => DAexec_time)
    
    return DAoutcome

end



@everywhere function RTNminus1proxy(NetworkData, MicroData, outB, DAoutcome)

    OPFNetworkData = prepareOPFNetworkData(NetworkData)

    Btype = OPFNetworkData[:Btype]
    Gbusid = OPFNetworkData[:Gbusid]
    Gtype = OPFNetworkData[:Gtype]
    RampDn = OPFNetworkData[:RampDn]
    RampUp = OPFNetworkData[:RampUp]
    fuelType = OPFNetworkData[:fuelType]
    Pmin = OPFNetworkData[:Pmin]
    Pmax = OPFNetworkData[:Pmax]
    Cost = OPFNetworkData[:Cost]
    Dbusid = OPFNetworkData[:Dbusid]
    VoLL = OPFNetworkData[:VoLL]
    FrBusid = OPFNetworkData[:FrBusid]
    ToBusid = OPFNetworkData[:ToBusid]
    x = OPFNetworkData[:x]
    fmax = OPFNetworkData[:fmax]
    STfmax = OPFNetworkData[:STfmax]
    BranchOutRate = OPFNetworkData[:BranchOutRate]
    N = OPFNetworkData[:N]
    G = OPFNetworkData[:G]
    D = OPFNetworkData[:D]
    L = OPFNetworkData[:L]
    G_n = OPFNetworkData[:G_n]
    D_n = OPFNetworkData[:D_n]
    Fr_n = OPFNetworkData[:Fr_n]
    To_n = OPFNetworkData[:To_n]
    hpd = OPFNetworkData[:hpd]
    dpw = OPFNetworkData[:dpw]
    wpy = OPFNetworkData[:wpy]
    hpy = hpd * dpw * wpy

    HCap = MicroData[:HCap]
    Pload = MicroData[:Pload]
    uDA = fillWithZeros(DAoutcome[:uDAsparse], G, hpd)
    PDA = fillWithZeros(DAoutcome[:PDAsparse], G, hpd)
    for g=1:G, h=1:hpd
        PDA[g,h] = (uDA[g,h] == 1)? min(max(PDA[g,h], Pmin[Gtype[g]]), Pmax[Gtype[g]]) : PDA[g,h]
    end

    #initialize array of branch availability from outage schedule
    aBranch = ones(Int, L)
    aBranch[outB] = 0

    #N - 1 contingencies
    C = L
    a = ones(Float64, L, C) - eye(L, C)
    prob = Array{Float64}(C)
    BranchOutRate_h = BranchOutRate / hpy
    for c=1:C
        prob[c] = 1.0
        for l=1:L
            ForcedOutProb = BranchOutRate_h[l] * exp(-BranchOutRate_h[l])
            dummy = (l == c)?ForcedOutProb:1-ForcedOutProb
            prob[c] *= dummy
        end
    end

    Dt = 60
    dt = 15

    RToutcome = Dict()
    RTcost_arr = []
    PRTdiff0_mat = Float64[]
    PRTdiffc_mat = [Float64[] for c=1:C]
    PRTShed0_mat = Float64[]
    PRTShedc_mat = [Float64[] for c=1:C]
    RTexec_time_arr = Float64[]

    for h=1:hpd

        tic()

        m = Model(solver = CplexSolver(CPX_PARAM_SIMDISPLAY=0, CPX_PARAM_THREADS=16, CPX_PARAM_EPGAP=0.05))

        @variable(m, Pup0[1:G] >= 0)
        @variable(m, Pdn0[1:G] >= 0)
        @variable(m, Pupc[1:G,1:C] >= 0)
        @variable(m, Pdnc[1:G,1:C] >= 0)
        @variable(m, f0[1:L])
        @variable(m, th0[1:N])
        @variable(m, STfc[1:L,1:C])
        @variable(m, STthc[1:N,1:C])
        @variable(m, fc[1:L,1:C])
        @variable(m, thc[1:N,1:C])
        @variable(m, PShed0[1:D] >= 0)
        @variable(m, PShedc[1:D,1:C] >= 0)

        for d=1:D
            setupperbound(PShed0[d], 0.0)
            for c=1:C
                setupperbound(PShedc[d,c], 0.0)
            end
        end

        @objective(m, Min, sum(Cost[Gtype[g]][end] * (Pup0[g] - Pdn0[g]) * 100 for g=1:G) )


        @constraint(m, Reference0, th0[1] == 0 )
        @constraint(m, ReferenceSTc[c=1:C], STthc[1,c] == 0 )
        @constraint(m, Referencec[c=1:C], thc[1,c] == 0 )

        #Pre-contingency state
        @constraint(m, Power0[n=1:N], 
            sum(PDA[g,h] + (Pup0[g] - Pdn0[g]) for g in G_n[n]) - 
            (sum(f0[l] for l in Fr_n[n]) - sum(f0[l] for l in To_n[n])) == 
            sum(Pload[d,h] - PShed0[d] for d in D_n[n])
            )
        @constraint(m, Flow0[l=1:L], f0[l] - aBranch[l] * (1/x[l]) * (th0[FrBusid[l]] - th0[ToBusid[l]]) == 0 )

        @constraint(m, Flow0PosLim[l=1:L], f0[l] <= fmax[l] )
        @constraint(m, Flow0NegLim[l=1:L], -f0[l] <= fmax[l] )

        @constraint(m, GenMin0[g=1:G], PDA[g,h] - Pdn0[g] >= uDA[g,h] * Pmin[Gtype[g]])
        @constraint(m, GenMax0[g=1:G; fuelType[Gtype[g]] != "Hydro"], PDA[g,h] + Pup0[g] <= uDA[g,h] * Pmax[Gtype[g]])
        @constraint(m, GenMax0Hydro[g=1:G; fuelType[Gtype[g]] == "Hydro"], PDA[g,h] + Pup0[g] <= uDA[g,h] * HCap[h])

        @constraint(m, RampDown0[g=1:G], Pdn0[g] <= 0.5 * Dt * RampDn[Gtype[g]] )
        @constraint(m, RampUp0[g=1:G], Pup0[g] <= 0.5 * Dt * RampUp[Gtype[g]] )

        #Short-term post-contingency state
        @constraint(m, STPowerc[n=1:N,c=1:C], 
            sum(PDA[g,h] + (Pup0[g] - Pdn0[g]) for g in G_n[n]) - 
            (sum(STfc[l,c] for l in Fr_n[n]) - sum(STfc[l,c] for l in To_n[n])) == 
            sum(Pload[d,h] - PShed0[d] - PShedc[d,c] for d in D_n[n])  
            )
        @constraint(m, STFlowc[l=1:L,c=1:C], STfc[l,c] - aBranch[l] * a[l,c] * (1/x[l]) * (STthc[FrBusid[l],c] - STthc[ToBusid[l],c]) == 0 )

        @constraint(m, STFlowcPosLim[l=1:L,c=1:C], STfc[l,c] <= STfmax[l])
        @constraint(m, STFlowcNegLim[l=1:L,c=1:C], -STfc[l,c] <= STfmax[l])

        #Post-contingency state (after the successful application of corrective control)
        @constraint(m, Powerc[n=1:N,c=1:C], 
            sum(PDA[g,h] + (Pup0[g] - Pdn0[g]) + (Pupc[g,c] - Pdnc[g,c]) for g in G_n[n]) - 
            (sum(fc[l,c] for l in Fr_n[n]) - sum(fc[l,c] for l in To_n[n])) == 
            sum(Pload[d,h] - PShed0[d] - PShedc[d,c] for d in D_n[n])
            )
        @constraint(m, Flowc[l=1:L,c=1:C], fc[l,c] - aBranch[l] * a[l,c] * (1/x[l]) * (thc[FrBusid[l],c] - thc[ToBusid[l],c]) == 0 )

        @constraint(m, FlowcPosLim[l=1:L,c=1:C], fc[l,c] <= fmax[l])
        @constraint(m, FlowcNegLim[l=1:L,c=1:C], -fc[l,c] <= fmax[l])

        @constraint(m, GenMinc[g=1:G,c=1:C], PDA[g,h] + Pup0[g] - Pdn0[g] - Pdnc[g,c] >= uDA[g,h] * Pmin[Gtype[g]])
        @constraint(m, GenMaxc[g=1:G,c=1:C; fuelType[Gtype[g]] != "Hydro"], PDA[g,h] + Pup0[g] - Pdn0[g] + Pupc[g,c] <= uDA[g,h] * Pmax[Gtype[g]])
        @constraint(m, GenMaxcHydro[g=1:G,c=1:C; fuelType[Gtype[g]] == "Hydro"], PDA[g,h] + Pup0[g] - Pdn0[g] + Pupc[g,c] <= uDA[g,h] * HCap[h])

        @constraint(m, RampDownc[g=1:G,c=1:C], Pdnc[g,c] <= dt * RampDn[Gtype[g]] )
        @constraint(m, RampUpc[g=1:G,c=1:C], Pupc[g,c] <= dt * RampUp[Gtype[g]] )

        TT = STDOUT # save original STDOUT stream
        redirect_stdout()
        status = solve(m; suppress_warnings=false)
        redirect_stdout(TT) # restore STDOUT

        if status == :Optimal
            append!(PRTdiff0_mat, getvalue(Pup0) - getvalue(Pdn0))
            append!(PRTShed0_mat, zeros(Float64, D))
            exp_corr_cost = 0.0
            for c=1:C
                corr_cost = 0.0
                for g=1:G
                    corr_cost += Cost[Gtype[g]][end] * (getvalue(Pupc)[g,c] - getvalue(Pdnc)[g,c]) * 100
                end
                exp_corr_cost += prob[c] * corr_cost
                append!(PRTdiffc_mat[c], [getvalue(Pupc)[g,c] - getvalue(Pdnc)[g,c] for g=1:G])
                append!(PRTShedc_mat[c], zeros(Float64, D))
            end
            push!(RTcost_arr, getobjectivevalue(m) + exp_corr_cost)
            push!(RTexec_time_arr, toq())
            continue
        end

        #activate load shedding
        for d=1:D
            setupperbound(PShed0[d], Pload[d,h])
            for c=1:C
                setupperbound(PShedc[d,c], Pload[d,h])
            end
        end

        @constraint(m, maxLoadShedding[d=1:D,c=1:C], Pload[d,h] - PShed0[d] - PShedc[d,c] >= 0)
        @objective(m, Min, (C+2) * sum(PShed0[d] for d=1:D) + sum(sum(PShedc[d,c] for d=1:D) for c=1:C)) 

        TT = STDOUT # save original STDOUT stream
        redirect_stdout()
        status = solve(m; suppress_warnings=false)
        redirect_stdout(TT) # restore STDOUT

        if status == :Optimal
            append!(PRTdiff0_mat, getvalue(Pup0) - getvalue(Pdn0))
            append!(PRTShed0_mat, getvalue(PShed0))
            prev_cost = 0.0
            for g=1:G
                prev_cost += Cost[Gtype[g]][end] * (getvalue(Pup0)[g] - getvalue(Pdn0)[g]) * 100
            end
            for d=1:D
                prev_cost += VoLL[d] * getvalue(PShed0)[d] * 100
            end
            exp_corr_cost = 0.0
            for c=1:C
                corr_cost = 0.0
                for g=1:G
                    corr_cost += Cost[Gtype[g]][end] * (getvalue(Pupc)[g,c] - getvalue(Pdnc)[g,c]) * 100
                end
                for d=1:D
                    corr_cost += VoLL[d] * getvalue(PShedc)[d,c] * 100
                end
                exp_corr_cost += prob[c] * corr_cost
                append!(PRTdiffc_mat[c], [getvalue(Pupc)[g,c] - getvalue(Pdnc)[g,c] for g=1:G])
                append!(PRTShedc_mat[c], [getvalue(PShedc)[d,c] for d=1:D])
            end
            push!(RTcost_arr, prev_cost + exp_corr_cost)
            push!(RTexec_time_arr, toq())
            continue
        end

        #shed all load
        append!(PRTdiff0_mat, zeros(Float64, G))
        append!(PRTShed0_mat, Pload[:,h])
        for c=1:C
            append!(PRTdiffc_mat[c], zeros(Float64, G))
            append!(PRTShedc_mat[c], zeros(Float64, D))
        end
        push!(RTcost_arr, 100 * sum(VoLL .* Pload[:,h]))
        push!(RTexec_time_arr, toq())

    end

    push!(RToutcome, :RTcost => RTcost_arr)
    push!(RToutcome, :PRTdiff0sparse => sparse(reshape(PRTdiff0_mat, G, hpd)))
    push!(RToutcome, :PRTdiffcsparsedict => Dict(zip(collect(1:C), [sparse(reshape(PRTdiffc_mat[c], G, hpd)) for c=1:C])))
    push!(RToutcome, :PRTShed0sparse => sparse(reshape(PRTShed0_mat, D, hpd)))
    push!(RToutcome, :PRTShedcsparsedict => Dict(zip(collect(1:C), [sparse(reshape(PRTShedc_mat[c], D, hpd)) for c=1:C])))
    push!(RToutcome, :RTexec_time => RTexec_time_arr)

    return RToutcome

end


@everywhere function assessmentProxy(NetworkData, MicroData, outB, DAoutcome, RToutcome, contingency_kwrd)

    OPFNetworkData = prepareOPFNetworkData(NetworkData)

    Btype = OPFNetworkData[:Btype]
    Vsp = OPFNetworkData[:Vsp]
    Bgs = OPFNetworkData[:Bgs]
    Bbs = OPFNetworkData[:Bbs]
    Gbusid = OPFNetworkData[:Gbusid]
    Gtype = OPFNetworkData[:Gtype]
    RampDn = OPFNetworkData[:RampDn]
    RampUp = OPFNetworkData[:RampUp]
    fuelType = OPFNetworkData[:fuelType]
    Pmin = OPFNetworkData[:Pmin]
    Pmax = OPFNetworkData[:Pmax]
    Qmin = OPFNetworkData[:Qmin]
    Qmax = OPFNetworkData[:Qmax]
    Cost = OPFNetworkData[:Cost]
    Dbusid = OPFNetworkData[:Dbusid]
    VoLL = OPFNetworkData[:VoLL]
    FrBusid = OPFNetworkData[:FrBusid]
    ToBusid = OPFNetworkData[:ToBusid]
    r = OPFNetworkData[:r]
    x = OPFNetworkData[:x]
    b = OPFNetworkData[:b]
    tr = OPFNetworkData[:tr]
    fmax = OPFNetworkData[:fmax]
    BranchOutRate = OPFNetworkData[:BranchOutRate]
    N = OPFNetworkData[:N]
    G = OPFNetworkData[:G]
    D = OPFNetworkData[:D]
    L = OPFNetworkData[:L]
    G_n = OPFNetworkData[:G_n]
    D_n = OPFNetworkData[:D_n]
    Fr_n = OPFNetworkData[:Fr_n]
    To_n = OPFNetworkData[:To_n]
    hpd = OPFNetworkData[:hpd]
    dpw = OPFNetworkData[:dpw]
    wpy = OPFNetworkData[:wpy]
    hpy = hpd * dpw * wpy
    BranchOutRate_h = BranchOutRate / hpy

    #compute admittance data
    for i in eachindex(tr)
        if tr[i] == 0.0 tr[i] = 1.0 end
    end
    tm = tr.^2
    gl = r ./ (r.^2 + x.^2)
    bl = -x ./ (r.^2 + x.^2)
    G_frfr = gl./tm
    G_frto = -gl./tr
    G_tofr = -gl./tr
    G_toto = gl
    B_frfr = (bl + 0.5*b)./tm
    B_frto = -bl./tr
    B_tofr = -bl./tr
    B_toto = bl + 0.5*b

    #micro-scenario data
    HCap = MicroData[:HCap]
    Pload = MicroData[:Pload]
    Qload = MicroData[:Qload]

    #DA and RT data
    uDA = fillWithZeros(DAoutcome[:uDAsparse], G, hpd)
    PDA = fillWithZeros(DAoutcome[:PDAsparse], G, hpd)
    PRTdiff0 = fillWithZeros(RToutcome[:PRTdiff0sparse], G, hpd)
    PRTdiffc = Dict()
    for (key,value) in RToutcome[:PRTdiffcsparsedict]
        push!(PRTdiffc, key => fillWithZeros(value, G, hpd))
    end
    PRTShed0 = fillWithZeros(RToutcome[:PRTShed0sparse], D, hpd)
    PRTShedc = Dict()
    for (key,value) in RToutcome[:PRTShedcsparsedict] 
        push!(PRTShedc, key => fillWithZeros(value, D, hpd))
    end

    #branch availability
    aBranch0 = ones(Int, L)
    aBranch0[outB] = 0

    #populate contingency list
    if contingency_kwrd == "N"
        contingency_list = [0]
    end
    if contingency_kwrd == "N-1"
        contingency_list = collect(1:L)
    end

    #indices of sync condensers
    syncs = find([Gtype[g] for g=1:G] .== "SynCond")

    #assessment
    assessment = zeros(Float64, hpd)
    assessment_time = zeros(Float64, hpd)

    for cont in contingency_list

        #compute contingency probability
        prob = 1.0
        for l=1:L
            ForcedOutProb = BranchOutRate_h[l] * exp(-BranchOutRate_h[l])
            dummy = (l in cont)?ForcedOutProb:1-ForcedOutProb
            prob *= dummy
        end

        #apply contingency
        aBranch = copy(aBranch0)
        for c in cont
            if c != 0 aBranch[c] = 0 end
        end        

        #apply preventive acctions
        ugen_spec = uDA
        Pgen_spec = PDA + PRTdiff0
        Pload_spec = Pload - PRTShed0

        #force sync condensors on
        ugen_spec[syncs] = 1

        #apply corrective actions if any
        if cont in keys(PRTdiffc) 
            Pgen_spec += PRTdiffc[cont] 
        end
        if cont in keys(PRTShedc)
            Pload_spec -= PRTShedc[cont] 
        end

        #force Pgen_spec within limits
        for g=1:G, h=1:hpd
            Pgen_spec[g,h] = (ugen_spec[g,h] == 1)? min(max(Pgen_spec[g,h], Pmin[Gtype[g]]), Pmax[Gtype[g]]) : Pgen_spec[g,h]
        end

        #check if actions are consistent, i.e. if the total generation is equal to the total load
        Pgen_tot = sum(Pgen_spec, 1)
        Pload_tot = sum(Pload_spec, 1)
        @assert maximum(abs.(Pgen_tot - Pload_tot)) < 5e3

        for h=1:1
#=
            println("Generator data [gen number|real power in pu]")
            for g=1:G
                if ugen_spec[g,h] == 0 continue end
                println("$g, $(ugen_spec[g,h]), $(Pgen_spec[g,h])")
            end

            println("Load data [load number|real pwr demand (pu)|reac. pwr demand (pu)]")
            for d=1:D
                println("$(Dbusid[d]), $(Pload_spec[d,h]), $(Qload[d,h])")
            end
=#

            tic()

            m = Model(solver = IpoptSolver(print_level=0))

            @variable(m, V[1:N])
            @variable(m, th[1:N])
            @variable(m, Pgen[1:G])
            @variable(m, Qgen[1:G])
            @variable(m, fp_fr[1:L])
            @variable(m, fp_to[1:L])
            @variable(m, fq_fr[1:L])
            @variable(m, fq_to[1:L])

            for n=1:N
                setvalue(V[n], 1.0)
                setvalue(th[n], 0.0)
            end

            @constraint(m, SWbusVoltage[n=1:N; Btype[n] == "SW"], V[n] == Vsp[n])
            @constraint(m, PVbusVoltageMax[n=1:N; Btype[n] == "PV"], V[n] <= Vsp[n] + 1e-5)
            @constraint(m, PVbusVoltageMin[n=1:N; Btype[n] == "PV"], V[n] >= Vsp[n] - 1e-5)
            @constraint(m, SWbusAngle[n=1:N; Btype[n] == "SW"], th[n] == 0)
            @constraint(m, PGenLosses[g=1:G; Btype[Gbusid[g]] != "SW"], Pgen[g] == Pgen_spec[g,h] )

            @NLconstraint(m, RealPowerBalance[n=1:N], 
                sum(Pgen[g] for g in G_n[n]) - (sum(fp_fr[l] for l in Fr_n[n]) + sum(fp_to[l] for l in To_n[n])) - Bgs[n] * V[n]^2 == 
                sum(Pload_spec[d,h] for d in D_n[n])
                ) 
            @NLconstraint(m, ReactivePowerBalance[n=1:N], 
                sum(Qgen[g] for g in G_n[n]) - (sum(fq_fr[l] for l in Fr_n[n]) + sum(fq_to[l] for l in To_n[n])) + Bbs[n] * V[n]^2 == 
                sum(Qload[d,h] for d in D_n[n])
                )

            @NLconstraint(m, RealFlowFromBus[l=1:L], fp_fr[l] - aBranch[l] * (G_frfr[l] * V[FrBusid[l]]^2 + (G_frto[l] * cos(th[FrBusid[l]] - th[ToBusid[l]]) + B_frto[l] * sin(th[FrBusid[l]] - th[ToBusid[l]])) * V[FrBusid[l]] * V[ToBusid[l]]) == 0 )
            @NLconstraint(m, RealFlowToBus[l=1:L], fp_to[l] - aBranch[l] * (G_toto[l] * V[ToBusid[l]]^2 + (G_tofr[l] * cos(th[ToBusid[l]] - th[FrBusid[l]]) + B_tofr[l] * sin(th[ToBusid[l]] - th[FrBusid[l]])) * V[ToBusid[l]] * V[FrBusid[l]]) == 0 )
            @NLconstraint(m, ReactiveFlowFromBus[l=1:L], fq_fr[l] - aBranch[l] * (-B_frfr[l] * V[FrBusid[l]]^2 + (G_frto[l] * sin(th[FrBusid[l]] - th[ToBusid[l]]) - B_frto[l] * cos(th[FrBusid[l]] - th[ToBusid[l]])) * V[FrBusid[l]] * V[ToBusid[l]]) == 0 )
            @NLconstraint(m, ReactiveFlowToBus[l=1:L], fq_to[l] - aBranch[l] * (-B_toto[l] * V[ToBusid[l]]^2 + (G_tofr[l] * sin(th[ToBusid[l]] - th[FrBusid[l]]) - B_tofr[l] * cos(th[ToBusid[l]] - th[FrBusid[l]])) * V[ToBusid[l]] * V[FrBusid[l]]) == 0 )


            @constraint(m, PQbusVoltageMax[n=1:N; Btype[n] == "PQ"], V[n] <= 1.1)
            @constraint(m, PQbusVoltageMin[n=1:N; Btype[n] == "PQ"], V[n] >= 0.9)
            @constraint(m, GenQmin[g=1:G; Btype[Gbusid[g]] != "SW"], Qgen[g] >= ugen_spec[g,h] * Qmin[Gtype[g]] )
            @constraint(m, GenQmax[g=1:G; Btype[Gbusid[g]] != "SW"], Qgen[g] <= ugen_spec[g,h] * Qmax[Gtype[g]] )
            @constraint(m, FlowFromPosLim[l=1:L], fp_fr[l]^2 + fq_fr[l]^2 <= fmax[l]^2 )
            @constraint(m, FlowToPosLim[l=1:L], fp_to[l]^2 + fq_to[l]^2 <= fmax[l]^2 )


            TT = STDOUT # save original STDOUT stream
            redirect_stdout()
            status = solve(m; suppress_warnings=false)
            redirect_stdout(TT) # restore STDOUT

            if status == :Optimal
                println("cont $(cont), hour $(h), no problem detected")
                for g=1:G
                    println(ugen_spec[g], " ", Qmin[Gtype[g]], " ", getvalue(Qgen[g]), " ", Qmax[Gtype[g]])
                end
                println(maximum(getvalue(V)), " ", minimum(getvalue(V)))
                assessment_time[h] += toq()
                continue
            end

            #activate load shedding
            println("activating load shedding")

            #compute participation factors
            totalPmax = sum([ugen_spec[g,h] * Pmax[Gtype[g]] for g=1:G])
            alpha = [ugen_spec[g,h] * Pmax[Gtype[g]] / totalPmax for g=1:G]

            m = Model(solver = IpoptSolver(print_level=0))

            @variable(m, V[1:N])
            @variable(m, th[1:N])
            @variable(m, Pgen[1:G])
            @variable(m, Qgen[1:G])
            @variable(m, fp_fr[1:L])
            @variable(m, fp_to[1:L])
            @variable(m, fq_fr[1:L])
            @variable(m, fq_to[1:L])
            @variable(m, PloadShed[1:D] >= 0)
            @variable(m, QloadShed[1:D])
            @variable(m, Ploss)

            for d=1:D
                setupperbound(PloadShed[d], Pload_spec[d,h])
            end

            for n=1:N
                setvalue(V[n], 1.0)
                setvalue(th[n], 0.0)
            end 

            @objective(m, Min, sum(PloadShed[d] for d=1:D) ) 

            @constraint(m, SWbusAngle[n=1:N; Btype[n] == "SW"], th[n] == 0)

            @constraint(m, PGenShed[g=1:G], Pgen[g] == Pgen_spec[g,h] - alpha[g] * (sum(PloadShed[d] for d=1:D) + Ploss) )

            @NLconstraint(m, RealPowerBalance[n=1:N], 
                sum(Pgen[g] for g in G_n[n]) - (sum(fp_fr[l] for l in Fr_n[n]) + sum(fp_to[l] for l in To_n[n])) - Bgs[n] * V[n]^2 == 
                sum(Pload_spec[d,h] - PloadShed[d] for d in D_n[n])
                ) 
            @NLconstraint(m, ReactivePowerBalance[n=1:N], 
                sum(Qgen[g] for g in G_n[n]) - (sum(fq_fr[l] for l in Fr_n[n]) + sum(fq_to[l] for l in To_n[n])) + Bbs[n] * V[n]^2 == 
                sum(Qload[d,h] - QloadShed[d] for d in D_n[n])
                )

            @NLconstraint(m, RealFlowFromBus[l=1:L], fp_fr[l] - aBranch[l] * (G_frfr[l] * V[FrBusid[l]]^2 + (G_frto[l] * cos(th[FrBusid[l]] - th[ToBusid[l]]) + B_frto[l] * sin(th[FrBusid[l]] - th[ToBusid[l]])) * V[FrBusid[l]] * V[ToBusid[l]]) == 0 )
            @NLconstraint(m, RealFlowToBus[l=1:L], fp_to[l] - aBranch[l] * (G_toto[l] * V[ToBusid[l]]^2 + (G_tofr[l] * cos(th[ToBusid[l]] - th[FrBusid[l]]) + B_tofr[l] * sin(th[ToBusid[l]] - th[FrBusid[l]])) * V[ToBusid[l]] * V[FrBusid[l]]) == 0 )
            @NLconstraint(m, ReactiveFlowFromBus[l=1:L], fq_fr[l] - aBranch[l] * (-B_frfr[l] * V[FrBusid[l]]^2 + (G_frto[l] * sin(th[FrBusid[l]] - th[ToBusid[l]]) - B_frto[l] * cos(th[FrBusid[l]] - th[ToBusid[l]])) * V[FrBusid[l]] * V[ToBusid[l]]) == 0 )
            @NLconstraint(m, ReactiveFlowToBus[l=1:L], fq_to[l] - aBranch[l] * (-B_toto[l] * V[ToBusid[l]]^2 + (G_tofr[l] * sin(th[ToBusid[l]] - th[FrBusid[l]]) - B_tofr[l] * cos(th[ToBusid[l]] - th[FrBusid[l]])) * V[ToBusid[l]] * V[FrBusid[l]]) == 0 )


            @constraint(m, VoltageMax[n=1:N], V[n] <= 1.05)
            @constraint(m, VoltageMin[n=1:N], V[n] >= 0.95)
            @constraint(m, GenPmin[g=1:G], Pgen[g] >= ugen_spec[g,h] * Pmin[Gtype[g]] )
            @constraint(m, GenPmax[g=1:G; fuelType[Gtype[g]] != "Hydro"], Pgen[g] <= ugen_spec[g,h] * Pmax[Gtype[g]] )
            @constraint(m, GenPmaxHydro[g=1:G; fuelType[Gtype[g]] == "Hydro"], Pgen[g] <= ugen_spec[g,h] * HCap[h] )
            @constraint(m, GenQmin[g=1:G], Qgen[g] >= ugen_spec[g,h] * Qmin[Gtype[g]] )
            @constraint(m, GenQmax[g=1:G], Qgen[g] <= ugen_spec[g,h] * Qmax[Gtype[g]] )
            @constraint(m, FlowFromPosLim[l=1:L], fp_fr[l]^2 + fq_fr[l]^2 <= fmax[l]^2 )
            @constraint(m, FlowToPosLim[l=1:L], fp_to[l]^2 + fq_to[l]^2 <= fmax[l]^2 )

            TT = STDOUT # save original STDOUT stream
            redirect_stdout()
            status = solve(m; suppress_warnings=false)
            redirect_stdout(TT) # restore STDOUT

            if status == :Optimal
                assessment[h] += prob * sum(VoLL .* (100 * getvalue(PloadShed)))
                assessment_time[h] += toq()
                println("cont $(cont), hour $(h), problem solved with minimal load-shedding, cost = $(prob * sum(VoLL .* (100 * getvalue(PloadShed))))")
                continue
            end

            #shed all load
            println("cont $(cont), hour $(h), problem solved with full load-shedding, cost = $(prob * sum(VoLL .* (100 * Pload_spec[:,h])))")
            assessment[h] += prob * sum(VoLL .* (100 * Pload_spec[:,h]))
            assessment_time[h] += toq()

        end

    end 


    return assessment, assessment_time

end

@everywhere function readOrComputeOutageOutcome(path, s, day, outB)

    sort!(outB)

    RToutcome = readRToutcome(path, s, day, outB)
    DAoutcome = readDAoutcome(path, s, day, outB)
    if nothing in values(RToutcome)
        NetworkData = readNetworkData(path)
        MicroData = readMicroData(path, s, day)
        if nothing in values(DAoutcome)
            DAoutcome = DANminus1proxy(NetworkData, MicroData, outB)
            writeDAoutcome(path, s, day, outB, DAoutcome)
        end
        RToutcome = RTNminus1proxy(NetworkData, MicroData, outB, DAoutcome)
        writeRToutcome(path, s, day, outB, RToutcome)
    else
        if nothing in values(DAoutcome)
            NetworkData = readNetworkData(path)
            MicroData = readMicroData(path, s, day)
            DAoutcome = DANminus1proxy(NetworkData, MicroData, outB)
            writeDAoutcome(path, s, day, outB, DAoutcome)
        end
    end

    return DAoutcome, RToutcome

end



@everywhere function readOrComputeOutageAssessmentOutcome(path, s, day, outB)

    sort!(outB)

    DAoutcome, RToutcome = readOrComputeOutageOutcome(path, s, day, outB)

    NassessOutcome = nothing #readAssessOutcome(path, s, day, outB, "N")
    if NassessOutcome == nothing
        NetworkData = readNetworkData(path)
        MicroData = readMicroData(path, s, day)
        Nminus1AssessOutcome = assessmentProxy(NetworkData, MicroData, outB, DAoutcome, RToutcome, "N")
        #writeAssessOutcome(path, s, day, outB, Nminus1AssessOutcome)
    end
#=
    Nminus1AssessOutcome = nothing #readAssessOutcome(path, s, day, outB, "N-1")
    if Nminus1AssessOutcome == nothing
        NetworkData = readNetworkData(path)
        MicroData = readMicroData(path, s, day)
        Nminus1AssessOutcome = assessmentProxy(NetworkData, MicroData, outB, DAoutcome, RToutcome, "N-1")
        #writeAssessOutcome(path, s, day, outB, Nminus1AssessOutcome)
    end
=#


    return RToutcome, DAoutcome #, NassessOutcome, Nminus1AssessOutcome, CMNminus2AssessOutcome

end
