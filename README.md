# Replication code for "Geopolitical Fragmentation and Trade"

This repository contains all the code that is required to replicate the baseline **estimation** and all **simulations** in the following article

Rodolfo G. Campos, Julia Estefania-Flores, Davide Furceri, Jacopo Timini,
"Geopolitical fragmentation and trade", *Journal of Comparative Economics*, 2023, Volume 51, Issue 4, December 2023, pp. 1289-1315, https://doi.org/10.1016/j.jce.2023.06.008.



Computational requirements
---------------------------

### Software Requirements

- **Stata**: the code was last run with Stata version 17. It has the following dependencies
  - `ftools`
  - `ppmlhdfe`
  - `reghdfe`
- **Julia**: the code was last run with Julia version 1.8. It has the following dependencies
  - `CSV`
  - `DataFrames`
  - `LinearAlgebra`
  - `Statistics`
- **Python**: the code was last run with Python version 3.10. It has the the following dependencies
  - `geopandas`
  - `pandas`
  - `numpy`
  - `matplotlib`


Data Availability
----------------------------

### Summary of Availability

- The database necessary for the estimation is in a zipped file called `geofragtrade.zip`. It is placed in the folder called `source_data`.
- All data necessary for the simulations are placed in the folder called `source_data`.
- The folder `map` contains the Natural Earth map, which is in the public domain. Visit the web site https://www.naturalearthdata.com/ for more details.

Further details
----------------------------

- Expected results for the simulations are saved in the folder `results`.

Instructions
---------------------------
1. To perform the estimation, unzip the file `geofragtrade.zip`, modify the path in the Stata dofile called `estimation.do` so that it points to the location of the database, and execute the dofile.
2. To perform the simulations, execute the Julia code `simulation_trade_fragmentation.jl`.
3. To produce the figures that appear in the paper, execute the Python codes `fig0.py`, `fig1.py`, `fig2.py`, `figA.py`. 
