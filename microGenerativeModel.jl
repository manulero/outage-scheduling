#!~/Downloads/julia-0.4.5/usr/bin/julia
@everywhere using JuMP
@everywhere using CPLEX
@everywhere using Distributions
@everywhere using OffsetArrays


@everywhere	include("./dataInputOutput.jl")

@everywhere function safesplit(data)

	return try split(data, " ") catch; data end

end


@everywhere function generateGenOutages(PQfactor, AnnualPeakLoadData, HourlyLoadFactor, GeneratorData, GenLookupData, weeks_per_year, days_per_week, hours_per_day, day0, day1)

	tic()

	days_per_year = weeks_per_year * days_per_week
	hours_per_week = hours_per_day * days_per_week 
	Gtype = GeneratorData["Type"]
	Pmax = Dict(zip(GenLookupData["Type"], [sum(float(safesplit(GenLookupData["Pout"][i]))) for i=1:length(GenLookupData["Pout"])]))
	GenMaintPeriods = Dict(zip(GenLookupData["Type"], GenLookupData["MaintPeriods"]))
	PAnnualPeakLoad = PQfactor * AnnualPeakLoadData["Pload"][1]
	PHourlyLoad = PAnnualPeakLoad * HourlyLoadFactor	
	outG = [Int64[] for day=1:days_per_year]
	G = length(Gtype)
	Gcap = [tuple(g, Pmax[Gtype[g]]) for g=1:G]
	orDeredGcap = sort(collect(Gcap), by=x->x[2], rev=true)
	VRLoad = [sum(PHourlyLoad[(week-1) * hours_per_week + 1:week * hours_per_week]) for week=1:weeks_per_year]

	#main loop
	for (g,cap) in orDeredGcap
		n = GenMaintPeriods[Gtype[g]]
		tempVRLoad = VRLoad
		max = maximum(tempVRLoad)
		for m = 1:n
			#pick the week with the lowest virtual residual load
			week = indmin(tempVRLoad)
			#schedule the maintenance
			startday = (week-1) * days_per_week + 1
			endday = startday + days_per_week - 1
			for day=startday:endday
				push!(outG[day], g)
			end
			#make sure the same week doesn't get picked again
			tempVRLoad[week] = tempVRLoad[week] + max
			#increase VRLoad by the capacity of g in the scheduled periods
			VRLoad[week] =  VRLoad[week] + cap
		end

	end

	println("Generator outages generated in $(toq()) seconds.")

	return OffsetArray(outG[day0:day1], day0:day1)

end


@everywhere function generateHydroCapacity(PQfactor, AnnualPeakLoadData, HourlyLoadFactor, HydroData, weeks_per_year, days_per_week, hours_per_day, day0, day1)

	tic()

	days_per_year = weeks_per_year * days_per_week
	hours_per_year = days_per_year * hours_per_day
	hours_per_quarter = round(Int, hours_per_year/4)
	HHourCap = zeros(Float64, hours_per_year)	

	PAnnualPeakLoad = PQfactor * AnnualPeakLoadData["Pload"][1]
	PHourlyLoad = PAnnualPeakLoad * HourlyLoadFactor
	HQuarterCap = HydroData["QuarterCap"]
	HQuarterEn = HydroData["QuarterEnergy"]

	#initialize first and last hour of the quarter
	first = 1
	last = hours_per_quarter
	for quarter = 1:4
		#create array of tuples (hour, load) for the quarter
		HourLoad = [tuple(h, PHourlyLoad[h]) for h=first:last]
		#sort the array by load in decreasing orDer
		orDHourLoad = sort(collect(HourLoad), by=x->x[2], rev=true)
		#compute capacity and available energy for the quarter
		quarterCap = HQuarterCap[quarter]
		quarterEn = HQuarterEn[quarter]
		for (h, load) in orDHourLoad
			#compute capacity to be assigned as the minimum between the capacity of the quarter and the energy available
			cap = min(quarterCap, quarterEn)
			#assign capacity
			HHourCap[h] = cap
			#reduce available energy
			quarterEn = quarterEn - cap
			if quarterEn <= 0 break end
		end
		#compute first and last hour of the next quarter
		first = first + hours_per_quarter
		last = last + hours_per_quarter

	end

	HCap = OffsetArray([HHourCap[(day - 1)*hours_per_day + 1: day*hours_per_day] for day=day0:day1], day0:day1)

	println("Hydro capacity generated in $(toq()) seconds.")

	return HCap

end


@everywhere function generateDemandForecast(PQfactor, AnnualPeakLoadData, HourlyLoadFactor, LoadData, hours_per_day, day0, day1)

	tic()

	PBusLoad = .01 * LoadData["P_PCT"]
	QBusLoad = .01 * LoadData["Q_PCT"]
	PAnnualPeakLoad = PQfactor * AnnualPeakLoadData["Pload"][1]
	QAnnualPeakLoad = PQfactor * AnnualPeakLoadData["Qload"][1]
	PHourlyLoad = PAnnualPeakLoad * HourlyLoadFactor
	QHourlyLoad = QAnnualPeakLoad * HourlyLoadFactor

	Ploadfc = OffsetArray([PBusLoad * transpose(PHourlyLoad[(day - 1)*hours_per_day + 1: day*hours_per_day]) for day=day0:day1], day0:day1)
	Qloadfc = OffsetArray([QBusLoad * transpose(QHourlyLoad[(day - 1)*hours_per_day + 1: day*hours_per_day]) for day=day0:day1], day0:day1)

	println("Demand forecast generated in $(toq()) seconds.")

	return Ploadfc, Qloadfc

end

@everywhere function generateDemandRealisation(LoadData, Ploadfc, Qloadfc, hours_per_day, day0, day1)

	tic()

	Pload = OffsetArray(Array{Float64}, day0:day1)
	Qload = OffsetArray(Array{Float64}, day0:day1)

	PBusLoad = .01 * LoadData["P_PCT"]
	QBusLoad = .01 * LoadData["Q_PCT"]
	D = length(PBusLoad)

	#Laurine's method
	sigma_glob = .03
	sigma_loc = 1.01 * sigma_glob
	Psum_of_squares = sum(PBusLoad.^2)
	Qsum_of_squares = sum(QBusLoad.^2)

	Palpha_d = Float16(sqrt((sigma_glob^2 - Psum_of_squares * sigma_loc^2) / (1 - Psum_of_squares)))
	Pbeta_d = Float16(sqrt((sigma_glob^2 - sigma_loc^2) / (Psum_of_squares - 1)))
	Qalpha_d = Float16(sqrt((sigma_glob^2 - Qsum_of_squares * sigma_loc^2) / (1 - Qsum_of_squares)))
	Qbeta_d = Float16(sqrt((sigma_glob^2 - sigma_loc^2) / (Qsum_of_squares - 1)))

	PN_alpha = Normal(0, Palpha_d)
	PN_beta = Normal(0, Pbeta_d)
	QN_alpha = Normal(0, Qalpha_d)
	QN_beta = Normal(0, Qbeta_d)

	for day = day0:day1
		Ptemp = []
		Qtemp = []
		for hour = 1:hours_per_day
			Peps_alpha = rand(PN_alpha)
			Peps_beta = rand(PN_beta, D)
			Qeps_alpha = rand(QN_alpha)
			Qeps_beta = rand(QN_beta, D)
			append!(Ptemp, (1 + Peps_alpha + Peps_beta) .* Ploadfc[day][:,hour])
			append!(Qtemp, (1 + Qeps_alpha + Qeps_beta) .* Qloadfc[day][:,hour])
		end
		Pload[day] = reshape(Ptemp, D, hours_per_day)
		Qload[day] = reshape(Qtemp, D, hours_per_day)
	end

	println("Demand realisation generated in $(toq()) seconds.")

	return Pload, Qload

end

@everywhere function generateMCOutcome(GeneratorData, GenLookupData, outG, Ploadfc, HCap, hpd, day0, day1)

	tic()

	uMCout = OffsetArray(Any, day0:day1)
	PMCout = OffsetArray(Any, day0:day1)
	vMCout = OffsetArray(Any, day0:day1)

	Gtype = GeneratorData["Type"]
	RampDn = Dict(zip(GenLookupData["Type"], GenLookupData["RampDn"]))
	RampUp = Dict(zip(GenLookupData["Type"], GenLookupData["RampUp"]))
	MinDnTime = Dict(zip(GenLookupData["Type"], GenLookupData["MinDnTime"]))
	MinUpTime = Dict(zip(GenLookupData["Type"], GenLookupData["MinUpTime"]))
	SUCost = Dict(zip(GenLookupData["Type"], GenLookupData["fuelPrice"] .* GenLookupData["SUCost"]))
	fuelType = Dict(zip(GenLookupData["Type"], GenLookupData["fuelType"]))
	Pout = Dict(zip(GenLookupData["Type"], [float(safesplit(GenLookupData["Pout"][i])) for i=1:length(GenLookupData["Pout"])]))
	Cost = Dict(zip(GenLookupData["Type"], GenLookupData["fuelPrice"] .* [float(safesplit(GenLookupData["Cost"][i])) for i=1:length(GenLookupData["Cost"])]))
	#set startup and shutdown ramp rates equal to the max generation
	RampSD = Dict(zip(GenLookupData["Type"], [sum(float(safesplit(GenLookupData["Pout"][i]))) for i=1:length(GenLookupData["Pout"])]))
	RampSU = RampSD	

	G = length(Gtype)
	D = size(Ploadfc[day0], 1)
	K = length(Pout[Gtype[1]])
	#initialize gen status and time of 'must keep status'
	uInit = zeros(Int, G)
	MinKeepTime = zeros(Int, G)

	Dt = 60


	for day=day0:day1

	    aGen = ones(Int, G)
	    aGen[outG[day]] = 0

		m = Model(solver = CplexSolver(CPX_PARAM_MIPDISPLAY=0, CPX_PARAM_EPGAP=0.001))

		@variable(m, PMC[1:G,1:hpd])
		@variable(m, PMCk[1:G,1:hpd,1:K] >= 0)
		@variable(m, uMC[1:G,1:hpd], Bin)
		@variable(m, 0 <= vMC[1:G,1:hpd] <= 1)
		@variable(m, 0 <= wMC[1:G,1:hpd] <= 1)

		@objective(m, Min, sum(sum(SUCost[Gtype[g]] * vMC[g,h] + sum(Cost[Gtype[g]][k] * PMCk[g,h,k] * 100 for k=1:K) for g=1:G) for h = 1:hpd) )

	    #Generator availability and cycling
		@constraint(m, keep[g=1:G, h = 1:MinKeepTime[g]], uMC[g,h] == uInit[g])
		@constraint(m, updown1[g=1:G], vMC[g,1] - wMC[g,1] == uMC[g,1] - uInit[g])
		@constraint(m, updown[g=1:G, h = 2:hpd], vMC[g,h] - wMC[g,h] == uMC[g,h] - uMC[g,h-1] )
		@constraint(m, MinUpTimeConstr[g=1:G, h = MinUpTime[Gtype[g]]:hpd], sum(vMC[g,s] for s=h-MinUpTime[Gtype[g]]+1:h) <= uMC[g,h] )
		@constraint(m, MinDnTimeConstr[g=1:G, h = MinDnTime[Gtype[g]]:hpd], sum(wMC[g,s] for s=h-MinDnTime[Gtype[g]]+1:h) <= 1 - uMC[g,h] )

		#Network balance and limits
		@constraint(m, NetworkBalance[h=1:hpd], sum(aGen[g] * PMC[g,h] for g=1:G) ==
			sum(Ploadfc[day][d,h] for d=1:D) )

		@constraint(m, MaxPLimk1[g=1:G, h=1:hpd], PMCk[g,h,1] == uMC[g,h] * Pout[Gtype[g]][1] )
		@constraint(m, MaxPLimk2to4[g=1:G, h=1:hpd, k=2:K], PMCk[g,h,k] <= uMC[g,h] * Pout[Gtype[g]][k] )
		@constraint(m, sumPk[g=1:G, h = 1:hpd], PMC[g,h] == sum(PMCk[g,h,k] for k=1:K))
		@constraint(m, HydroCap[g=1:G, h = 1:hpd; fuelType[Gtype[g]] == "Hydro"], PMCk[g,h,K] <= uMC[g,h] * HCap[day][h] )

		@constraint(m, RampDownConstr[g=1:G, h=2:hpd], PMC[g,h-1] - PMC[g,h] <= Dt * RampDn[Gtype[g]] * uMC[g,h] + RampSD[Gtype[g]] * wMC[g,h])
		@constraint(m, RampUpConstr[g=1:G, h=2:hpd], PMC[g,h] - PMC[g,h-1] <= Dt * RampUp[Gtype[g]] * uMC[g,h-1] + RampSU[Gtype[g]] * vMC[g,h])

		status = solve(m)

		uMCstar = round.(Int, getvalue(uMC))
		vMCstar = round.(Int, getvalue(vMC))
		PMCstar = getvalue(PMC)

		uMCout[day] = sparse(uMCstar)
		PMCout[day] = sparse(PMCstar)
		vMCout[day] = sparse(vMCstar)

	end

	println("MC outcome generated in $(toq()) seconds.")

	return uMCout, PMCout, vMCout

end

@everywhere function computeHourlyLoadFactor(WeeklyPeakLoadData, DailyPeakLoadData, HourlyPeakLoadData)

	WeeklyPeakLoad = .01 * WeeklyPeakLoadData["PCTLoad"]
	DailyPeakLoad = .01 * DailyPeakLoadData["PCTLoad"]
	HourlyPeakLoad = .01 * HourlyPeakLoadData["PCTLoadWinterWeekday"]

	weeks_per_year = length(WeeklyPeakLoad)
	days_per_week = length(DailyPeakLoad)
	hours_per_day = length(HourlyPeakLoad)

	HourlyLoadFactor = []

	for week = 1:weeks_per_year
		for day = 1:days_per_week
			for hour = 1:hours_per_day
        		day_type = day <= 5 ? "Weekday" : "Weekend"
				if week <= 8 || week >= 44
		            season = "Winter"
		        elseif week >= 18 && week <= 30
		            season = "Summer"
		        else
		            season = "Spring"
		        end
		        HourlyPeakLoad = .01 * HourlyPeakLoadData["PCTLoad$(season)$(day_type)"]
		        push!(HourlyLoadFactor, HourlyPeakLoad[hour] * DailyPeakLoad[day] * WeeklyPeakLoad[week])
		    end
		end
    end

    return HourlyLoadFactor, weeks_per_year, days_per_week, hours_per_day

end

@everywhere function genMicros(s, day0, day1, path)


	println("\nGenerating micro-scenarios\n")
	println("Days $(day0) to $(day1)")
	srand(s)

	BusData = readData("$(path)/BusData.csv")
	LoadData = readData("$(path)/LoadData.csv")
	GeneratorData = readData("$(path)/GeneratorData.csv")
	GenLookupData = readData("$(path)/GenLookupData.csv")
	BranchData = readData("$(path)/BranchData.csv")
	HydroData = readData("$(path)/HydroData.csv")
	AnnualPeakLoadData = readData("$(path)/AnnualPeakLoadData.csv")
	WeeklyPeakLoadData = readData("$(path)/WeeklyPeakLoadData.csv")
	DailyPeakLoadData = readData("$(path)/DailyPeakLoadData.csv")
	HourlyPeakLoadData = readData("$(path)/HourlyPeakLoadData.csv")

	HourlyLoadFactor, weeks_per_year, days_per_week, hours_per_day = computeHourlyLoadFactor(WeeklyPeakLoadData, DailyPeakLoadData, HourlyPeakLoadData)

	U = Uniform(0.9, 1.1)
	PQfactor = rand(U)

	outG = generateGenOutages(PQfactor, AnnualPeakLoadData, HourlyLoadFactor, GeneratorData, GenLookupData, weeks_per_year, days_per_week, hours_per_day, day0, day1)
	HCap = generateHydroCapacity(PQfactor, AnnualPeakLoadData, HourlyLoadFactor, HydroData, weeks_per_year, days_per_week, hours_per_day, day0, day1)
	Ploadfc, Qloadfc = generateDemandForecast(PQfactor, AnnualPeakLoadData, HourlyLoadFactor, LoadData, hours_per_day, day0, day1)
	Pload, Qload = generateDemandRealisation(LoadData, Ploadfc, Qloadfc, hours_per_day, day0, day1)
	uMC, PMC, vMC = generateMCOutcome(GeneratorData, GenLookupData, outG, Ploadfc, HCap, hours_per_day, day0, day1)
	
	writeMicroData(path, s, day0, day1, outG, HCap, Ploadfc, Qloadfc, Pload, Qload, uMC, PMC, vMC)

end

#=





@everywhere function generateFcdBranchOutages(BranchData, weeks_per_year, days_per_week, hours_per_day, day0, day1)

	tic()

	days_per_year = days_per_week * weeks_per_year
	hours_per_year = hours_per_day * days_per_year 	
	BranchOutRate = BranchData["BranchOutRate"] / hours_per_year
	BranchOutDur = round.(Int, BranchData["BranchOutDur"])
	L = length(BranchOutRate)

	BranchOutages = [[] for l=1:L]

	for l=1:L
		E = Exponential(1 / BranchOutRate[l])
		clock = 0
		while(true)
			fail = clock + ceil(Int, rand(E))
			if fail > hours_per_year break end
			rep = fail + BranchOutDur[l] - 1
			push!(BranchOutages[l], (fail, rep))
			clock = rep
		end
	end

	FcdOutB = [zeros(Int, L, hours_per_day) for day=1:days_per_year]

	for l=1:L
		for (fail, rep) in BranchOutages[l]
			for hour = fail:min(rep, hours_per_year)
				day, hour_of_day = ceil(Int, hour / hours_per_day), 1 + ((hour - 1) % hours_per_day)
				FcdOutB[day][l, hour_of_day] = 1
			end
		end
	end

	println("Branch forced outages generated in $(toq()) seconds.")

	return OffsetArray([sparse(FcdOutB[day]) for day=day0:day1], day0:day1)

end

### CONTINGENCIES
#GENERATOR FORCED OutAGES (IN HOURS)

tic()
#first create outage array of arrays 
GenForceOut = Array(Any, G)

for g=1:G
	#Compute hourly outage rate (outages per hour)
	hGenOutRate = GenOutRate[Gtype[g]] / (hpd * 365)
	#Compute hourly outage probability
	hGenOutProb = hGenOutRate * exp(-hGenOutRate)
	#Compute outage expected duration (in hours)
	ExpGenOutDur = GenOutDuration[Gtype[g]]
	#create array that will be filled with the tuples 'outage time, duration' in hours
	GenForceOut[g] = Tuple{Int,Int}[]
	#start sweeping all hours in the year
	hour = 1
	while hour <= hpy
		day = Int(ceil(hour/24))
		week = Int(ceil(day/7))
		#check if generator already in scheduled maintenance for current week
		if week in GenSchOut[g]
			hour = hour + hpw	#go to next week
			continue
		end
		#simulate forced outage
		if rand(Uniform()) <= hGenOutProb
			#simulate outage duration
			GenOutDur = round(Int, rand(Exponential(ExpGenOutDur)))
			#add tuple
			push!(GenForceOut[g], (hour, GenOutDur))
			hour = hour + GenOutDur
			continue
		end
		hour = hour + 1	
	end
end
push!(T, toq())
println("Task ", ARGS[2], ": Generator forced outages took ", T[end])


#BRANCH FORCED OutAGES (IN HOURS)

tic()
#first create outage array of arrays 
BranchForceOut = Array(Any, L)
#Compute hourly outage rate (outages per hour)
hBranchOutRate = BranchOutRate / hpy

for l=1:L
	#Compute hourly outage probability
	hBranchOutProb = hBranchOutRate[l] * exp(-hBranchOutRate[l])
	#create array that will be filled with the tuples 'outage time, duration' in hours
	BranchForceOut[l] = Tuple{Int, Int}[]
	#start sweeping all hours in the year
	hour = 1
	while hour <= hpy
		day = Int(ceil(hour/24))
		#assume scheduled maintenance is not decided yet		
		#simulate forced outage
		if rand(Uniform()) <= hBranchOutProb
			#simulate outage duration
			BranchOutDur = round(Int, rand(Normal(BranchOutDuration[l], 1)))
			#add tuple
			push!(BranchForceOut[l], (hour, BranchOutDur))
			hour = hour + BranchOutDur
			continue
		end
		hour = hour + 1	

	end
end
push!(T, toq())
println("Task ", ARGS[2], ": Branch forced outages took ", T[end])
=#