# Brayton-Cycle-Variable-cp
Thermodynamic modelling of a gas turbine Brayton cycle using MATLAB to compare constant and temperature-dependent specific heat effects on efficiency and performance.
## Numerical Analysis of Brayton Cycle Performance Using Temperature-Dependent Specific Heat
This project investigates how temperature-dependent specific heat affects the performance prediction of a gas turbine Brayton cycle. Traditional thermodynamic analysis assumes constant specific heat (cp), however at high turbine inlet temperatures this assumption becomes inaccurate. A MATLAB numerical model was developed to compare constant cp and variable cp formulations based on NASA thermodynamic data.
## Objectives
- Develop a MATLAB numerical model of the Brayton cycle using constant and variable specific heat (cp) formulations
- Use entropy balance equations to calculate compressor and turbine exit temperatures for the variable cp model
- Evaluate performance parameters: net work output, thermal efficiency, back work ratio and optimum pressure ratio
- Perform a parametric study by varying compressor pressure ratio and turbine inlet temperature
- Quantify the error introduced by the constant cp assumption at high turbine inlet temperatures
- Compare performance trends of constant and variable cp models and explain the physical reasons for the differences
## Methodology
- Ideal Brayton cycle analysis
- NASA polynomial thermodynamic property relations
- Numerical integration of cp(T)
- Iterative root-finding solution for state temperatures
- Parametric study for pressure ratio 5–50 and TIT 1200–1800 K
## Key Results
- Constant cp model over-predicts thermal efficiency at high temperatures
- Optimum pressure ratio shifts when variable cp is used
- Back work ratio increases when real gas behaviour is considered
- Error increases significantly as turbine inlet temperature rises
## Software Used
- MATLAB
- Numerical integration
- Root-finding algorithms
## Engineering Interpretation
The increase in specific heat at high temperature increases compressor work and reduces turbine expansion effectiveness. Therefore simplified Brayton cycle analysis gives unrealistic performance predictions for modern gas turbines. Variable specific heat modelling provides more physically realistic cycle behaviour.
## Future Work
Future work will include component efficiencies, pressure losses, intercooling and regeneration to better represent real gas turbine engines.
## Project Structure
- /code – MATLAB implementation of the Brayton cycle model
- /plots – performance comparison graphs
- README.md – project description and results summary
## How to Run
1. Open MATLAB
2. Download the repository files
3. Open the file brayton_cycle_variable_cp.m
4. Run the script
5. The program calculates compressor and turbine states and generates performance plots
---
## Sample Results
### Work Output Comparison
![Work Plot](plots/.png)
