# OWF O&M OPTIMIZATION
Thesis repository containing OWF O&M model created by Jonathan Hjorthøj-Nielsen and Johan Østergaard Knarreborg. 

<div align="center">
  <img width="400" src="https://github.com/JohanKnarreborg/OWF_OM_OPTIMIZATION/blob/setup/figs/maintenance.gif">
  <br>
  <em>Corrective Maintenance (blue) and Scheduled Maintenance (white) for a yearly simulation for Kriegers Flak💨</em>
</div>

## Project description 
One of the main cost drivers of offshore wind farms (OWF) is the Operations and Maintenance (O&M) activities conducted throughout the lifetime of a wind farm. As seen in (Röckmann, Lagerveld, and Stavenuiter 2017) the cost of O&M related activities can account for up to 25%-30% of the total lifetime costs, and therefore crucial to minimize.

This study aims to reduce the O&M related costs by building a decision making model that can make informed decisions based on forecasts of power prices as well as power production. This will allow an OWF operator to reduce the opportunity costs of shutting down their wind turbines to conduct maintenance and reducing the lifetime cost of the OWF. 


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
The output from each model can be explored using the output analysis notebook. Here example plots for an output analysis is available. 

#### Model output  
The model outputs a number of files containing information about all model decisions. The model outputs are: 

| **File**     | **Type** | **Size**                                 | **Description**                                                                                   |
|--------------|----------|------------------------------------------|---------------------------------------------------------------------------------------------------|
| **r_sm.npy** | Bin      | [quarters, turbine]                      | Whether each turbine has scheduled maintenance performed in each quarter.                         |
| **h_sm.npy** | Int      | [quarters, turbine]                      | How many quarters of scheduled maintenance is remaining for each turbine.                         |
| **r_cm.npy** | Bin      | [quarters, turbine]                      | Whether each turbine has corrective maintenance performed in each quarter.                        |
| **h_cm.npy** | Int      | [quarters, turbine]                      | How many quarters of corrective maintenance is remaining for each turbine.                        |
| **Z.npy**    | Bin      | [CTV, quarters, nodes]                   | Whether each CTV is at node n in quarter t.                                                       |
| **MC.npy**   | Int      | [days]                                   | Daily to total cost of fuel.                                                                      |
| **X.npy**    | Bin      | [quarters, turbine]                      | Whether each turbine is producing in each quarter. Rather use r_sm and h_cm to determine if turbine is running. This will be zero when no wind is forecasted or zero prices are forecasted. |

## Input data
TBD...

##  Model architecture
The OR model is a mixed-integer linear programming model that optimizes the O&M activities for an OWF. The model is based on a rolling horizon approach where the model is solved for each day in the simulation. For each day a part of the OR model is redefined to fit the current day. The basic model running structure is shown below: 

<div align="center">
  <img width="600" src="https://github.com/JohanKnarreborg/OWF_OM_OPTIMIZATION/blob/main/figs/flow_chart.png">
  <br>
  <em>Model architecture for long-term simulations. </em>
</div>

##  Code structure 
The main OR model is developed in julia. All other code for this project has been developed in python 3.12. The code is structured as follows:
    
```
README.md               # This file
requirements.txt        # Requirements file for python environment
model                   # Main model folder
|   model.jl            # Julia file containing the model
|   output              # Output folder from the model
|   |   ...
|   data                # Folder for storing main data for model 
|   |   layout          # Folder containing all data related to OWF layout
|   |   |   ...
|   |   maintenance     # Folder containing the failure data for turbines
|   |   |   ...
|   |   prices          # Folder containing all price forecasts
|   |   |   ...
|   |   sets            # Folder containing the sets used by the model
|   |   |   ... 
|   |   weather         # Folder containing weather forecast and actual weather 
|   |   |   ... 
|   |   jl_function.jl  # Julia file containing the functions used by model
|   data_generation     # Folder for data generation
|   |   failures        # Generation file for failures 
|   |   |   ...
|   |   price_forecasts # Price forecast files 
|   |   |   ... 
|   |   weather_forecast# Weather forecasts generating the sunthetic forecast
|   |   |   ...
output_analysis         # Folder for output analysis
|   output_analysis.ipynb
```
