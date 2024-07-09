using JuMP, Gurobi, HiGHS 
using CSV
using NPZ
using ArgParse
using DataFrames
using Dates
include("data/julia_function.jl")


#define timing variables to time model definition and solving time
global timing_model_definition = 0 
global timing_model_solving = 0    
timing_total = @elapsed begin
timing_model_definition_inter = @elapsed begin


#Adding run arguments 
s = ArgParseSettings()
@add_arg_table! s begin
    "--p_cost"
    help = "Fixed postponement cost"
    arg_type = Int
    default = 75

    "--days"
    help = "Number of days to run"
    arg_type = Int
    default = 365

    "--w"
    help = "Number of turbines"
    arg_type = Int
    default = 36

    "--e_price"
    help = "Electricity price for a fixed price forecast. (Production optimized)"
    arg_type = Int
    default = 0

    "--day_start"
    help = "Start day"
    arg_type = Int
    default = 1

    "--price_quantile"
    help = "Quantile, e.g., 80. Not used if 0. "
    arg_type = Int
    default = 0

    "--comment"
    help = "Comment to add to output file"
    arg_type = String
    default = ""

    "--h_sm_per_turbine"
    help = "Quarters of SM per turbine"
    arg_type = Int
    default = 160

    "--heuristic_turbines"
    help = "Number of turbines used in the matheuristic"
    arg_type = Int
    default = 5

    "--forecast"
    help = "Forecast method"
    arg_type = String
    default = "EWMA"

    "--ctvs"
    help = "Number of CTVs"
    arg_type = Int
    default = 1

    "--crews"
    help = "Number of crews"
    arg_type = Int
    default = 4

    "--naive_method_quarters"
    help = "Number of quarters the NAIVE method should do each day. If 0, the NAIVE method is not used."
    arg_type = Int
    default = 0

    "naive_method_startday"
    help = "Day to start the naive method"
    arg_type = Int
    default = 5

    "--year"
    help = "Year"
    arg_type = Int
    default = 2021

    "--co2tax"
    help = "CO2 tax if 0 normal fuel cost is used."
    arg_type = Int
    default = 0

    "--solver"
    help = "Either Gurobi or HiGHS solver"
    arg_type = String
    default = "HiGHS"
end

parsed_args = parse_args(ARGS, s)

# Assign parsed arguments to variables
solver = parsed_args["solver"]
PP_cost_input = parsed_args["p_cost"]
CTV = parsed_args["ctvs"]
days = parsed_args["days"]
W = parsed_args["w"]
N = W + CTV 
CREW = parsed_args["crews"]
CTV_CREW_CAP = CREW+1
argH_sm = zeros(W)
argH_sm .= parsed_args["h_sm_per_turbine"]  
argH_cm = zeros(W)
DAY_start_initial = parsed_args["day_start"]
year = parsed_args["year"]
price_quantile = parsed_args["price_quantile"]
naive_method_quarters = parsed_args["naive_method_quarters"]
naive_method_startday = parsed_args["naive_method_startday"]
heuristic_turbines = parsed_args["heuristic_turbines"]
forecast = parsed_args["forecast"]
e_price = parsed_args["e_price"] 
CO2_tax = parsed_args["co2tax"]
comment = parsed_args["comment"] #add comment to output file 

CM_cost = 1000
T = 96 #number of quarters to look ahead. 96 if planning horizon is one day. 192 if two days.
path_pre = "" #add a path prefix if needed
DAY_end = 24*4 #ending time of day
FC = 0.25 * (1+CO2_tax/100) #0.25 is the price of fuel
T_hours = convert(Int64,T/4) ## number of hours used 
Time_period_lookout = 6 # nr of periods that the model looks into the future , max 12, min 3 (since CTV needs to be able to reach first wind turbine.)
WHL = 1.8 #wave height limit
WSL = 20 #wind speed limit
PP_cost = PP_cost_input #cost of postponement
M_cm = 12*4*14*10 #BigM for CM 
M_wave_wind_limit = 100 #BigM for wave and wind limit
H_sm = zeros(Int64, T, W) #scheduled maintenance hours for each wind turbine at each time period
H_cm = zeros(Int64, T, W) #corrective maintenance hours for each wind turbine at each time period


######## Optimizer attributes #########################
if solver == "Gurobi"
    m = Model(Gurobi.Optimizer,add_bridges = false) 
    set_optimizer_attributes(m, "OutputFlag" => 1)
    set_optimizer_attribute(m, "MIPGap", 0.00)
    set_optimizer_attribute(m, "Seed", 42*2) 
    #set_optimizer_attribute(m,"NodefileStart", 0.5) #0-inf
    #set_optimizer_attribute(m,"Presolve", 2) #0,1,2
    #set_optimizer_attribute(m, "Threads", 16) 
    #set_optimizer_attribute(m, "Method", 3)
    #set_optimizer_attribute(m, "BarOrder", 0)
else
    m = Model(HiGHS.Optimizer)
end

#########         SETS        #########################
#O as type int
O = CSV.read(path_pre*"data/sets/O.csv", DataFrame, header=0) #night quarters set 

#########       DATA         #########################
#weather data
if year == 2021 
    wind_forecast_mat = npzread(path_pre*"data/weather/2021_wind_forecast.npy")
    wh_forecast_mat =CSV.read(path_pre*"data/weather/2021_data.csv", DataFrame, header=1)
    if price_quantile != 0
        prices_mat =CSV.read(path_pre*"data/weather/2021_data.csv", DataFrame, header=1)
    end

elseif year == 2022 
    wind_forecast_mat = npzread(path_pre*"data/weather/2022_wind_forecast.npy")
    wh_forecast_mat =CSV.read(path_pre*"data/weather/2022_data.csv", DataFrame, header=1)
    if price_quantile != 0
        prices_mat =CSV.read(path_pre*"data/weather/2022_data.csv", DataFrame, header=1)
    end

else 
    println("Error: No forecast")
end



#price forecast 
if forecast == "EWMA"
    if year == 2021
        p_forecast = CSV.read(path_pre*"data/prices/EWMA_2021.csv", DataFrame, header=1) 
    elseif year == 2022
        p_forecast = CSV.read(path_pre*"data/prices/EWMA_2022.csv", DataFrame, header=1) 
    end
    p_forecast_mat = p_forecast


elseif forecast == "ANN"
    if year == 2021
        p_forecast = CSV.read(path_pre*"data/prices/XGBoost_ANN_2021.csv", DataFrame, header=1) 
    elseif year == 2022
        p_forecast = CSV.read(path_pre*"data/prices/XGBoost_ANN_2022.csv", DataFrame, header=1)
    end
    p_forecast_mat = p_forecast


elseif forecast == "ARIMAX"
    if year == 2021
        p_forecast = CSV.read(path_pre*"data/prices/ARIMAX_2021.csv", DataFrame, header=1)
    elseif year == 2022
        p_forecast = CSV.read(path_pre*"data/prices/ARIMAX_2022.csv", DataFrame, header=1)
    end
    p_forecast_mat = p_forecast

elseif forecast == "PERFECT"
    if year == 2021
        p_forecast = CSV.read(path_pre*"data/prices/PERFECT_2021.csv", DataFrame, header=1)
    elseif year == 2022
        p_forecast = CSV.read(path_pre*"data/prices/PERFECT_2022.csv", DataFrame, header=1)
    end
    p_forecast_mat = p_forecast

else 
    #return error 
    println("Error: No forecast")
    
end

#if production based create a fixed price forecast
if e_price != 0
    p_forecast_mat = e_price*ones(days,T_hours)
else
    p_forecast_mat = p_forecast_mat
end


#maintenance hours 
if year == 2021
    H_sm_input = CSV.read(path_pre*"data/maintenance/H_sm.csv", DataFrame)[1:365,:]
    H_cm_input = CSV.read(path_pre*"data/maintenance/H_cm.csv", DataFrame)[1:365,:]
elseif year == 2022
    H_sm_input = CSV.read(path_pre*"data/maintenance/H_sm.csv", DataFrame)[363:728,:]
    H_cm_input = CSV.read(path_pre*"data/maintenance/H_cm.csv", DataFrame)[363:728,:]
end


#########       Distance Matrix     #####################
#create travel matrix 
D_time = CSV.read(path_pre*"data/layout/dist_matrix_time.csv", DataFrame, header=0)
D_meter = CSV.read(path_pre*"data/layout/dist_matrix_meter.csv", DataFrame, header=0)
D_meter = D_meter[2:end,:]
D_time = D_time[2:end,:]
#create a custom distance matrix based on the number of turbines and CTVs
D = zeros(N,N)
D[1:W,1:W] .= D_time[1:W,1:W]
for ctv in 1:CTV
    for w in 1:W
        D[W+ctv,w] = D_time[73,w]
        D[w,W+ctv] = D_time[w,73]
    end
end

#create consumption matrix
con_matrix = npzread(path_pre*"data/layout/con_matrix.npy")
con_matrix = con_matrix[:,:,1:Time_period_lookout]
con = zeros(N,N,Time_period_lookout) 
con[1:W,1:W,:] .= con_matrix[1:W,1:W,:]
for ctv in 1:CTV
    for w in 1:W
        con[W+ctv,w,:] = con_matrix[73,w,:]
        con[w,W+ctv,:] = con_matrix[w,73,:]
    end
end
#round con matrix 
con = round.(con,digits = 0)



end #end of timing model definition inter
global timing_model_definition += timing_model_definition_inter


for day in 1:days 
    global timing_model_definition_inter = @elapsed begin 
    global PP_cost, T_delta, year, price_quantile, naive_method_quarters, heuristic_turbines, forecast, e_price, CM_cost,argH_cm,argH_sm, DAY_start_initial,D_meter,H_cm,H_sm,wind_forecast_extended

    DAY_start = DAY_start_initial + day -1

    print("Day: ")
    println(DAY_start)
    
    T_delta = (DAY_start-1)*24+1
    H_sm = zeros(Int64, T, W) #scheduled maintenance hours for each wind turbine at each time period
    H_cm = zeros(Int64, T, W) #maintenance hours for each wind turbine at each time period
    

    #set price forecast 
    global p_forecast = Vector(p_forecast_mat[DAY_start,1:T_hours])
    p_forecast_extended = repeat(p_forecast,inner = 4) 

    #set weather forecast
    wind_forecast = wind_forecast_mat[DAY_start,:]
    wind_forecast_extended = repeat(wind_forecast,inner = 4)

    #power forecast 
    power_forecast = wind_to_power(wind_forecast, 8.4e6, 12, 3, 25, 3.1415*(167/2)^2, 0.37)/1e6
    global power_forecast_extended = repeat(power_forecast,inner = 4) 
  
    #if NAIVE method then no power forecast
    if naive_method_quarters != 0 
        power_forecast_extended = ones(T)
    end

    #wh forecast
    wh_forecast = wh_forecast_mat[T_delta:T_delta+T_hours-1,7] 
    wh_forecast_extended = repeat(wh_forecast,inner = 4)

    #if price quantile is set then calculate the postponement cost. Can only be done after 45 days (1080 hours) if using past 6 weeks prices
    if DAY_start > 45 && price_quantile != 0 && e_price == 0
            prices = prices_mat[T_delta-1056:T_delta,2] 
            global PP_cost = quantile(prices,price_quantile/100)
    end

    

    #on first run the non-day specific model is defined. 
    if DAY_start == DAY_start_initial
        
        ############## update data      ##################### 
        H_sm[1,:] = argH_sm #update remaining hours of schedule maintenance for each wind turbine at start time
        H_sm[1,:] = H_sm[1,:] + Vector(H_sm_input[DAY_start,1:W]) #add scheduled maintenance hours

        H_cm[1,:] = argH_cm #update remaining hours of unscheduled maintenance for each wind turbine at start time
        H_cm[1,:] = H_cm[1,:] + Vector(H_cm_input[DAY_start,1:W]) #add unscheduled maintenance hours


        ############### VARIABLES       #####################
        
        println("variables start")
        @variable(m, x[1:T,1:W], Bin)

        @variable(m, MC[1:T] >=0) 
        @variable(m, PP >=0)

        @variable(m,r_cm[1:T,1:W], Bin) 
        @variable(m,r_sm[1:T,1:W], Bin) 
        @variable(m,h_cm[1:T,1:W] >= 0,Int) 
        @variable(m,h_sm[1:T,1:W] >= 0,Int) 

        @variable(m, Z[1:CTV,1:T,1:N], Bin) #location of CTV i at time t, pos 73 is habor
        @variable(m, crew[1:T,1:CTV] >= 0, Int) #number of crews on each CTV at time period t

        @variable(m, MS_cm[1:T,1:W], Bin)
        @variable(m, ME_cm[1:T,1:W], Bin)
        @variable(m, MS_sm[1:T,1:W], Bin)
        @variable(m, ME_sm[1:T,1:W], Bin)
        @variable(m, MS[1:T,1:W], Bin) 
        @variable(m, ME[1:T,1:W], Bin)


        println("arcs start")
        @variable(m, A[s=1:N, t=1:T, r=1:N, tt=t:t+Time_period_lookout], Bin) #arc from node i at time t to node j at time t' 
        println("arcs done")


        ########## Objective function    ##################
        #pre-compute revenue for each quarter
        rev= [(power_forecast_extended[t] * p_forecast_extended[t])/4  for t in 1:T]
        @objective(m, Max, sum(x[t,w]*rev[t] for t in 1:T, w in 1:W) - sum(MC[t] for t = 1:T)  - PP )


        ######## CONSTRAINTS         ########################
        #maintenance constraints
        @constraint(m, [t=1:T],  MC[t] == sum(A[s,t,r,tt]*con[s,r,tt-t] for s=1:N,r=1:N,tt=t+1:t+Time_period_lookout)*FC  ) 

        @constraint(m, [t=1:T, w=1:W], h_cm[t,w] <= (1-x[t,w])*M_cm)
        @constraint(m, [t=1:T-1, w=1:W], h_cm[t+1,w] == h_cm[t,w]-r_cm[t,w])
        @constraint(m, [t=1:T-1, w=1:W], h_sm[t+1,w] == h_sm[t,w]-r_sm[t,w])

        @constraint(m, [t=1:T, w=1:W], (1-r_cm[t,w]) >= x[t,w])
        @constraint(m, [t=1:T, w=1:W], (1-r_sm[t,w]) >= x[t,w])

        @constraint(m, [t=2:T, w=1:W], MS_cm[t,w] >= r_cm[t,w]-r_cm[t-1,w])
        @constraint(m, [t=1:T-1, w=1:W], ME_cm[t,w] >= r_cm[t,w]-r_cm[t+1,w]) 
        @constraint(m, [t=2:T, w=1:W], MS_sm[t,w] >= r_sm[t,w]-r_sm[t-1,w])
        @constraint(m, [t=1:T-1, w=1:W], ME_sm[t,w] >= r_sm[t,w]-r_sm[t+1,w])

        @constraint(m, [t=1:T, w=1:W], 1 >= r_sm[t,w]+r_cm[t,w])

        @constraint(m, [t=1:T, w=1:W], 2*MS[t,w] >= MS_cm[t,w]+MS_sm[t,w])  #2 is big M 
        @constraint(m, [t=1:T, w=1:W], 2*ME[t,w] >= ME_cm[t,w]+ME_sm[t,w])  #2 is big M 


        #nighttime/off hours
        for t in O[:,1][1:convert(Int64,T_hours/2)]
            t = convert(Int64,t)
            @constraint(m, [q=1:4,w=1:W],r_cm[q + 4*(t-1),w] == 0) #no maintenance at night
            @constraint(m, [q=1:4,w=1:W],r_sm[q + 4*(t-1),w] == 0) #no maintenance at night
            @constraint(m, [q=1:4,ctv=1:CTV],Z[ctv,q + 4*(t-1),N-(CTV-ctv)]  == 1) #each CTV in each respective habor at night
        end 

        

        println("arc constraints start")
        @constraint(m, [s=1:N,t=1:T,r = 1:N], A[s,t,r,t] == 0) #no jumping to same point in time 
        
        @constraint(m, [s=1:N,t=Time_period_lookout+1:T-Time_period_lookout], sum(A[r,tt,s,t] for r = 1:N, tt = t-Time_period_lookout:t-1) == sum(A[s,t,r,tt] for r = 1:N, tt = t+1:t+Time_period_lookout))#arcs in = arcs out 
        for i = 1:Time_period_lookout 
            @constraint(m, [s=1:N,t=[i+1]], sum(A[r,tt,s,t] for r = 1:N, tt = 1:t-1) == sum(A[s,t,r,tt] for r = 1:N, tt = t+1:t+Time_period_lookout)) #arcs in = arcs out 
        end

        @constraint(m, [s=1:N,t=1:T-Time_period_lookout,r=1:N,tt=t+1:t+Time_period_lookout], A[s,t,r,tt] <= max(0,(tt-t) - D[s,r]) ) #only allow arcs that are possible due to travel matrix   

        @constraint(m, [t=[1]], sum(A[s,t,r,tt] for s = 1:N,r=1:N,tt=t+1:t+Time_period_lookout) == CTV)  #only CTV number of arcs out of time period one

        @constraint(m, [t=1:T-Time_period_lookout,s = 1:N], sum(A[s,t,r,tt] for r=1:N,tt=t+1:t+Time_period_lookout) <= 1)  #only one arc out of a node at a time
        for i = 1:Time_period_lookout 
            @constraint(m, [t=[T-i],s = 1:N], sum(A[s,t,r,tt] for r=1:N,tt=t+1:t+i) <= 1)   #only one arc out of a node at a time
        end

        @constraint(m, [s=1:N,t=1:T-Time_period_lookout,tt=t+2:t+Time_period_lookout], A[s,t,s,tt] == 0 )  #prevent CTVs traveling to same node as it is in, in more than one timestep 
        println("arc constraints done")


        #arc to ctv placement
        @constraint(m, [r=1:N,tt=Time_period_lookout+1:T], sum(A[s,t,r,tt] for s=1:N,t = tt-Time_period_lookout:tt-1)  == sum(Z[ctv,tt,r] for ctv = 1:CTV)) #If arc enters node then there must be a CTV
        for i = 2:Time_period_lookout 
            @constraint(m, [r=1:N,tt=[i]], sum(A[s,t,r,tt] for s=1:N,t = 1:tt-1)  == sum(Z[ctv,tt,r] for ctv = 1:CTV))
        end

        #follow each CTV 
        @constraint(m,[t=1:T-Time_period_lookout,s=1:N,r=1:N,tt=t+1:t+Time_period_lookout,ctv=1:CTV], A[s,t,r,tt] +Z[ctv,t,s] <= Z[ctv,tt,r] +1 )

        #CTV placement to Maintenance 
        @constraint(m, [t=1:T, w=1:W], sum(Z[ctv,t,w] for ctv in 1:CTV) >= MS[t,w]) #IF maintenance ends or starts one CTV needs to be at location
        @constraint(m, [t=1:T, w=1:W], sum(Z[ctv,t,w] for ctv in 1:CTV) >= ME[t,w]) #IF maintenance ends or starts one CTV needs to be at location

        #Crew placement 
        @constraint(m, [t=[1]], sum(crew[t,ctv] for ctv in 1:CTV) <= CREW) #Total number of crews
        @constraint(m, [t=1:T, ctv=1:CTV], crew[t,ctv] <= CTV_CREW_CAP) #Max number of crews on each CTV   
        @constraint(m, [t=2:T, ctv=1:CTV,w=1:W], Z[ctv,t,w] +  MS[t,w] - 1 <= (crew[t-1,ctv]-crew[t,ctv])   ) #Crew dropoff
        @constraint(m, [t=1:T-1, ctv=1:CTV,w=1:W], Z[ctv,t,w] +  ME[t,w] -1 <= (crew[t+1,ctv]-crew[t,ctv]) ) #Crew pickup
        @constraint(m, [t=1:T], sum(r_sm[t,w]+r_cm[t,w] for w in 1:W) + sum(crew[t,ctv] for ctv in 1:CTV) == CREW) #Total number of restricted in each time period
        println("crew constraints done")

        @constraint(m, [t=1:T, w=1:W], sum(Z[ctv,t,w] for ctv in 1:CTV) <= 1) #CTV only at one location at a time 
        @constraint(m, [ctv=1:CTV], sum(Z[ctv,t,W+ctv_other] for t in 1:T, ctv_other in 1:CTV if ctv_other != ctv) == 0) # this constraint prevents CTVs from switching harbours
        #for ctv in 1:CTV-1
        #    @constraint(m, sum(Z[ctv,t,W] for t in 1:T,w=1:W) >= sum(Z[ctv+1,t,W] for t in 1:T,w=1:W)) #Use CTVs in order
        #end


        println("CTV constraints done")

        ################performance constraints (not strictly neccessary)####################
        for o in O[:,1][1:convert(Int64,T_hours/2)] 
            o = convert(Int64,o)
            for q in 1:4
                t = q + 4*(o-1)
                @constraint(m, [s=1:N,r = 1:W, tt=t+1:t+Time_period_lookout],A[s,t,r,tt] == 0) #no arcs outside harbor at night
                @constraint(m, [s=1:W,r = 1:N, tt=t+1:t+Time_period_lookout],A[s,t,r,tt] == 0) #no arcs outside harbor at night

                #no ME and MS at night 
                @constraint(m, [w=1:W], ME[t,w] == 0) 
                @constraint(m, [w=1:W], MS[t,w] == 0)
                @constraint(m, [w=1:W], MS_cm[t,w] == 0)
                @constraint(m, [w=1:W], ME_cm[t,w] == 0)
                @constraint(m, [w=1:W], MS_sm[t,w] == 0)
                @constraint(m, [w=1:W], ME_sm[t,w] == 0)


            end
        end 

        #define upper limits
        @constraint(m, sum(MC) <= T*CTV*150*FC) 
        @constraint(m,[t =1:T,w=1:W], h_cm[t,w] <= 12*4*31*10) 
        @constraint(m,[t = 1:T,w =1:W], h_sm[t,w] <= 12*4*31)
        @constraint(m, [t = 1:T,ctv = 1:CTV],crew[t,ctv] <= CREW)
 
        ######################### END of initial setup ########################
    end 
    




    ####################### Deletion of constraints ############################
    if DAY_start >= DAY_start_initial + 1 
        global m,argH_cm,argH_sm, H_sm,H_cm,x, MC, PP, r_cm, r_sm, h_cm, h_sm, Z, crew, MS_cm, ME_cm, MS_sm, ME_sm, MS, ME, A, NAIVE_1, NAIVE_2, NAIVE_3,WHL_1, WHL_2, WHL_3, WSL_1, WSL_2, WSL_3, HEU_1, HEU_2, HEU_3, HEU_con, HSCH_1, HCAT_1, PP_con,naive_method_startday


        ############## update data      ##################### 
        H_sm[1,:] = argH_sm #update remaining hours of schedule maintenance for each wind turbine at start time
        H_sm[1,:] = H_sm[1,:] + Vector(H_sm_input[DAY_start,1:W]) #add scheduled maintenance hours

        H_cm[1,:] = argH_cm #update remaining hours of unscheduled maintenance for each wind turbine at start time
        H_cm[1,:] = H_cm[1,:] + Vector(H_cm_input[DAY_start,1:W]) #add unscheduled maintenance hours


        #new obj  
        rev= [(power_forecast_extended[t] * p_forecast_extended[t])/4  for t in 1:T]
        @objective(m, Max, sum(x[t,w]*rev[t] for t in 1:T, w in 1:W) - sum(MC[t] for t = 1:T)  - PP )


        #see this link on deletion times https://discourse.julialang.org/t/jump-delete-is-slower-than-building-a-new-model/92597/6. Constraints updated between runs should not be sparse arrays, since deleting these is very slow.
        delete(m, vec(WHL_1))
        delete(m, vec(WSL_1))
        delete(m, vec(HEU_1))
        delete(m, vec(HEU_2))
        delete(m, vec(HEU_3))
        delete(m, vec(HSCH_1))
        delete(m, vec(HCAT_1))
        delete.(m, PP_con)
        

        if (naive_method_quarters != 0)
            delete.(m, NAIVE_1)
        end
    end
    





    ##### Define day-specific constraints ####### 
    WHL_1 = @constraint(m,[t=1:T,ctv=1:CTV],wh_forecast_extended[t,1] <= WHL+ M_wave_wind_limit*(Z[ctv,t,N-(CTV-ctv)])) #in habour if WHL
    WSL_1 = @constraint(m,[t=1:T,ctv=1:CTV],wind_forecast_extended[t,1] <= WSL+ M_wave_wind_limit*(Z[ctv,t,N-(CTV-ctv)])) #in habour if WSL
    HCAT_1 = @constraint(m, [ w=1:W], h_cm[1,w] == H_cm[1,w])
    HSCH_1 = @constraint(m, [ w=1:W], h_sm[1,w] == H_sm[1,w])
    PP_con = @constraint(m, PP == sum(h_sm[T,w]*PP_cost/4 for w in 1:W) + CM_cost*sum(h_cm[T,:])) 


    ########################  HEURISTICS  definition ############################
    if heuristic_turbines == 0
        heuristic_turbines = CREW*2
    end

    turbines_with_maintenance = Int[]
    if sum(H_cm) > 0 #if CM exists then find heuristic_turbines closest turbines to the CM turbine and add to [turbines_with_maintenance]
        println("Cat maintenance needed")
        #find arg max of H_cm
        cat_index = argmax(H_cm[1,:]) 
        print("cat index: ")
        println(cat_index)
        D_meter_dist = Vector(D_meter[cat_index,1:W])
        D_meter_args = sortperm(D_meter_dist)
        for w in D_meter_args
            if length(turbines_with_maintenance) >= heuristic_turbines
                break
            end
            if sum(H_sm[:,w]) > 0
                push!(turbines_with_maintenance,w)
            end
        end


    else #if no CM exists add the closest turbines to harbour
        println("No cat maintenance needed")
        for w in 1:W
            #if length of turbines_with_maintenance is 5 then break 
            if length(turbines_with_maintenance) >= heuristic_turbines
                break
            end
            if sum(H_sm[:,w]) > 0
                push!(turbines_with_maintenance,w)
            end
        end 
    end
    println("turbines with maintenance: ")
    println(turbines_with_maintenance)


    turbine_list = zeros(Int, W)
    # Update the new list based on the existence of elements in the short list. This method is used due to the constraint deletion time of sparse arrays. 
    for w in 1:W
        if w in turbines_with_maintenance || sum(H_cm[:,w]) > 0
            turbine_list[w] = 1
        end
    end

    HEU_1 = @constraint(m, [t=1:T,ctv=1:CTV,w = 1:W ], Z[ctv,t,w] - turbine_list[w]  <= 0) 
    HEU_2 = @constraint(m, [t=1:T,w=1:W],MS[t,w] - turbine_list[w]  <= 0)
    HEU_3 = @constraint(m, [t=1:T,w=1:W],ME[t,w] - turbine_list[w]  <= 0)
    
    ########## NAIVE METHOD ############################
    if (naive_method_quarters != 0) 
        
        if (DAY_start < naive_method_startday) || (sum(wh_forecast_extended[32:84] .>= WHL) >= 0.01)  #if WHL or WSL is enforced or Naive start day not reached then no SM. 
            NAIVE_1 = @constraint(m, sum(r_sm[t,w] for t in 1:96, w in 1:W) == 0)

        else  
            #naive should do maximum of naive_method_quarters and sum of SM left for turbines included in heuristic_turbines
            if (sum(H_sm[1,w] for w in turbines_with_maintenance;init=0) <= naive_method_quarters)
                NAIVE_1 = @constraint(m, sum(r_sm[t,w] for t in 1:96, w in 1:W) == sum(H_sm[1,w] for w in turbines_with_maintenance;init=0))
            else 
                NAIVE_1 = @constraint(m, sum(r_sm[t,w] for t in 1:96, w in 1:W) == naive_method_quarters)
            
            end
        end

    end


    
    println(H_sm[1,:])
    println(H_cm[1,:])


    #end definition timing 
    end
    global timing_model_definition += timing_model_definition_inter

    ##################################################################
    println("optimizing")
    timing_model_solving_inter = @elapsed begin 
    optimize!(m)
    end
    global timing_model_solving += timing_model_solving_inter

    #define the amount of maintenance left to be used in the next iteration
    argH_sm = round.(JuMP.value.(h_sm)[DAY_end,:],digits = 0)
    argH_cm = round.(JuMP.value.(h_cm)[DAY_end,:],digits = 0)

    # export daily results. 
    global H_sm_ = round.(JuMP.value.(h_sm)[1:96,:],digits = 0)
    global h_cm_ = round.(JuMP.value.(h_cm)[1:96,:],digits = 0)
    global r_sm_ = round.(JuMP.value.(r_sm)[1:96,:],digits = 0) 
    global r_cm_ = round.(JuMP.value.(r_cm)[1:96,:],digits = 0) 
    global Z_ = round.(JuMP.value.(Z)[:,1:96,:],digits = 0) 
    global X_ = round.(JuMP.value.(x)[1:96,],digits = 0) 
    global MC_ = JuMP.value.(MC[1:96,])
    global MC_cost = sum(MC_[t] for t = 1:24*4) 
    

    #save daily results
    if DAY_start == DAY_start_initial
        global H_sm_ls = H_sm_
        global h_cm_ls = h_cm_
        global r_sm_ls = r_sm_
        global r_cm_ls = r_cm_
        global Z_ls = Z_
        global X_ls = X_
        global MC_ls = MC_cost
    else
        H_sm_ls = cat(H_sm_ls,H_sm_,dims=1)
        h_cm_ls = cat(h_cm_ls,h_cm_,dims=1) 
        r_sm_ls = cat(r_sm_ls,r_sm_,dims=1)
        r_cm_ls = cat(r_cm_ls,r_cm_,dims=1)
        Z_ls = cat(Z_ls,Z_,dims=2)  
        X_ls = cat(X_ls,X_,dims=1)
        MC_ls = cat(MC_ls,MC_cost,dims=1)
    end

end 
end 

println("Total time:")
println(timing_total)


output_folder = path_pre*"output/output_"*string(days)*"_"*string(W)*"_"*string(PP_cost_input)*"_"*string(e_price)*"_"*string(DAY_start_initial)*"_"*forecast*"_"*"_"*string(CREW)*"_"*string(CTV)*"_"*string(naive_method_quarters)*"_"*string(year)*"_"*string(price_quantile)*"_"*string(CO2_tax) *"_"*comment
#check if exists 
if !isdir(output_folder)
    mkdir(output_folder)
end
#save data as npy 
npzwrite(output_folder*"/h_sm.npy",H_sm_ls) 
npzwrite(output_folder*"/h_cm.npy",h_cm_ls)
npzwrite(output_folder*"/r_sm.npy",r_sm_ls)
npzwrite(output_folder*"/r_cm.npy",r_cm_ls)
npzwrite(output_folder*"/Z.npy",Z_ls)
npzwrite(output_folder*"/X.npy",X_ls)
npzwrite(output_folder*"/MC.npy",MC_ls)

#save yaml file 
yaml_file = output_folder*"/output.yaml"
open(yaml_file, "w") do io
    println(io, "Total time: $timing_total")
    println(io, "Model definition time: $timing_model_definition")
    println(io, "Model solving time: $timing_model_solving")
    println(io, "Days: $days")
    println(io, "Turbines: $W")
    println(io, "PP_cost: $PP_cost_input")
    println(io, "e_price: $e_price")
    println(io, "Fuel price: $FC")
    println(io, "Day_start: $DAY_start_initial")
    println(io, "Forecast: $forecast")
    println(io, "Crews: $CREW")
    println(io, "CTVs: $CTV")
    println(io, "naive_method_quarters: $naive_method_quarters")
    println(io, "naive_method_startday: $naive_method_startday")
    println(io, "Year: $year")
    println(io, "Price_quantile: $price_quantile")
    println(io, "Comment: $comment")
    println(io, "Timestamp: $(Dates.now())")
end








