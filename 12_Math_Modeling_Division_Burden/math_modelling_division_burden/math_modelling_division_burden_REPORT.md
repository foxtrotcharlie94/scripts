# math_modelling_division_burden.R

**Script:** `math_modelling_division_burden.R`

## What it does
Implements a stochastic branching-process simulation of HSC divisions (self-renewal, asymmetric/differentiation, death, fraction dividing per timestep) and compares simulated mean division burden and pool-size trajectories against closed-form analytical formulas; sweeps FRAC and A2 parameters and shows a burden-distribution histogram at the final timestep.

## Inputs
No external data files — pure parameter-driven simulation.

## Outputs
8 individual PNG/PDF plots plus a combined 8-panel HSC_panel.pdf/png.

## Conceptual aim / target
Provide a theoretical model of how division burden accumulates in an HSC population under different self-renewal/differentiation/death regimes, to interpret the empirical HB/LB clone phenotypes.

## Note
The only purely simulation/theory script in the set; byte-identical copy exists at Downloads root.
