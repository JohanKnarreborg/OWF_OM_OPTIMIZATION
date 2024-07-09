# OWF O&M OPTIMIZATION
Thesis repository containing OWF O&M model created by Jonathan Hjorthøj-Nielsen and Johan Knarreborg. 

<div align="center">
  <img width="400" src="https://github.com/JohanKnarreborg/OWF_OM_OPTIMIZATION/blob/setup/figs/maintenance.gif">
  <br>
  <em>Corrective Maintenance (blue) and Scheduled Maintenance (white) for a yearly simulation 🤯</em>
</div>

## Project description 
One of the main cost drivers of offshore wind farms (OWF) is the Operations and Maintenance (O\&M) activities conducted throughout the lifetime of a wind farm. As seen in \cite{OM_costs} the cost of O\&M related activities can account for up to 25\%-30\% of the total lifetime costs. 

Since the start of the energy crisis at the end of 2021 there has been an increase in electricity prices as well as volatility (\cite{energy_crisis}). This should further incentivize power producers to maximize their sale of power at price peaks and place maintenance activities in periods of low production and prices. By reducing the O\&M related costs from operating an OWF the business case of investing and owning OWF should be improved. This should allow an easier transition to a more sustainable future with renewable electricity sources. It should also provide a more certain investment opportunity as the costs from O\&M activities will be more predictable.

This study aims to reduce the O\&M related costs by building a decision making model that can make informed decisions based on forecasts of power prices as well as power production. This will allow an OWF owner to reduce the opportunity costs of shutting down their wind turbines to conduct maintenance. 


## 🚀 Quick start 🚀
The model can easlisy be run from the command line. First navigate to the model folder. Make sure you have julia installed and define your julia path. Most likely something like: '/Applications/Julia-1.8.app/Contents/Resources/julia/bin/julia'

Call the model.jl file with the applicable model parameters passed as arguments. If you do not have access to a Gurobi license, 'HiGHS' can be passed as the chosen solver. Here is an example model simulation using 12 turbines: 

'/Applications/Julia-1.8.app/Contents/Resources/julia/bin/julia' model.jl --p_cost=75 --days=31 --w=12 --e_price=0 --day_start=1  --comment=comment --h_sm_per_turbine=160 --heuristic_turbines=6 --forecast=EWMA --ctvs=1 --crews=2  --year=2021 --solver=Gurobi 

The output will be saved in the output folder in a subfolder named after the parameters used. 

## Model parameters 
--p_cost: (Int, default: 75)
Fixed postponement cost.

--days: (Int, default: 365)
Number of days to run the simulation.

--w: (Int, default: 36)
Number of turbines.

--e_price: (Int, default: 0)
Electricity price for a fixed price forecast. (Production optimized)

--day_start: (Int, default: 1)
Start day.

--price_quantile: (Int, default: 0)
Quantile, e.g., 80. If 0 then not used. 

--comment: (String, default: "")
Comment to add to output file.

--h_sm_per_turbine: (Int, default: 160)
Quarters of SM per turbine.

--heuristic_turbines: (Int, default: 5)
Number of turbines used in the matheuristic.

--forecast: (String, default: "EWMA")
Forecast method.

--ctvs: (Int, default: 1)
Number of CTVs.

--crews: (Int, default: 4)
Number of crews.

--naive_method_quarters: (Int, default: 0)
Number of quarters the NAIVE method should do each day. If 0, the NAIVE method is not used.

--naive_method_startday: (Int, default: 5)
Day to start the naive method.

--year: (Int, default: 2021)
Year.

--co2tax: (Int, default: 0)
CO2 tax; if 0, normal fuel cost is used.

--solver: (String, default: "HiGHS")
Either Gurobi or HiGHS solver.

## Output analysis 

The output from each model can be explored using the output analysis notebook.
