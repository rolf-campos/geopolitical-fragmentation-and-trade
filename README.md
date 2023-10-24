# Replication code for the simulations in "Geopolitical Fragmentation and Trade"

This repository contains all the code that is required to replicate the **simulations** in the following paper:

Rodolfo G. Campos, Julia Estefania-Flores, Davide Furceri, Jacopo Timini,
"Geopolitical fragmentation and trade", *Journal of Comparative Economics*, 2023, https://doi.org/10.1016/j.jce.2023.06.008.



Computational requirements
---------------------------

### Software Requirements

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

- All data necessary for the simulations are placed in the folder called `source_data`.
- The folder `map` contains the Natural Earth map, which is in the public domain. Visit the web site https://www.naturalearthdata.com/ for more details.

Further details
----------------------------

- Expected results are saved in the folder `results`.

Instructions
---------------------------
1. To perform the simulations, execute the Julia code `simulation_trade_fragmentation.jl`.
2. To produce the figures that appear in the paper, execute the Python codes `fig0.py`, `fig1.py`, `fig2.py`, `figA.py`. 
