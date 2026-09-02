# LPTTValidation

A Basilisk C solver Lptt.c for fully developed, axisymmetric pipe flow of a (L-PTT) viscoelastic fluid with a Newtonian solvent contribution.

The solver uses the log-conformation formulation (lptt-log-conform.h) to remain
numerically stable at high Weissenberg number, and is validated here against the
analytical velocity profiles of Cruz et al.[1] over a sweep of solvent
viscosity ratios, BETA = μ_s / μ_₀.

## Problem setup

- 2-D axisymmetric, periodic pipe flow driven by a constant body force.
- The relaxation time λ is set from the Weissenberg number, WI = λ·U_N/R,
  using the Newtonian mean velocity reference U_N = a_x·R²/(8·μ0) = 1/8.
- Domain size 1, grid resolution 2^7, no-slip at the wall (u.t[top] = dirichlet(0)).
- epsPTT = 0.25 sets the PTT nonlinearity parameter; trA = 3 initializes the
  conformation tensor.

## Files

 This solver includes lptt-log-conform.h (the log-conformation L-PTT header) .

 Make sure that header is present alongside Lptt.c.

## Compilation & execution

```bash
qcc -O2 -Wall -fopenmp Lptt.c -o run -lm

export OMP_NUM_THREADS=12     

> adjust to the machine being used

./run

```

The run advances to t = 70 (steady state) and writes the velocity profile to
profile_beta_BETA.dat as two columns: radial position y and normalized
axial velocity u.x / U_N.

To reproduce a different $β$ case (and the corresponding plot in Plots/), edit
BETA in Lptt.c to the desired value and recompile/rerun. The
Weissenberg number $Wi$ can be swept the same way.

## Validation plots

Plots/ contains the steady-state velocity profile for each BETA tested against
the Cruz et al. [1] analytical solution:

- validation_beta_0.01.png — BETA = 0.01
- validation_beta_0.1.png — BETA = 0.1
- validation_beta_0.5.png — BETA = 0.5
- validation_beta_0.8.png — BETA = 0.8
- validation_beta_1.png — BETA = 1 
- validation_all_beta.png — all cases overlaid

## Reference

[1] D. O. A. Cruz, F. T. Pinho, and P. J. Oliveira, "Analytical solutions for fully
developed laminar flow of some viscoelastic liquids with a Newtonian solvent
contribution," *Journal of Non-Newtonian Fluid Mechanics*, 132, 28–35, 2005.
