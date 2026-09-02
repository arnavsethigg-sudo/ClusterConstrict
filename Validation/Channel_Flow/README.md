# Droplet_Migration

A Basilisk C solver DropChannel.c for the axisymmetric migration of a
viscoelastic droplet through a straight cylindrical channel, using an
embedded boundary for the wall. The setup validates against the Newtonian
droplet migration and deformation results of Nath et al. (2017) [1].

The droplet's dispersed phase can carry a viscoelastic (L-PTT) contribution
via lptt-log-conform-embed.h, the embedded-boundary variant of the
log-conformation L-PTT header, which skips solid/cut cells belonging to the
embedded wall in every stress update.

## Problem setup

- 2-D axisymmetric channel of length 20, radius RIN = 1, resolved on a
  1<<9 base grid with adaptive refinement (adapt_wavelet on cs, f, u.x, u.y).
- The wall is embedded : wall_radius(x) currently returns a
  constant RIN, giving a straight channel; a curved/constricted geometry can
  be substituted by redefining wall_radius and the commented-out
  RTHROAT / XSTART / XEND parameters.
- A droplet of radius DROPLET_RADIUS = 0.8 is initialised at x = X0 = -8,
  upstream in the channel, and convects downstream under an imposed parabolic
  inlet velocity profile.
- Continuous and dispersed phases are neutrally buoyant (RHO_RATIO = 1) in
  this configuration.
- Control parameters, all non-dimensional:
  - RE — continuous-phase Reynolds number
  - CA — capillary number, sets surface tension from the continuous-phase
    viscosity
  - WI — Weissenberg number, sets the droplet-phase relaxation time LAM
  - MU_RATIO — ratio of total droplet viscosity to continuous-phase viscosity
  - BETA — solvent fraction of the droplet-phase viscosity
- epsPTT = sets the PTT nonlinearity parameter; trA is initialised to 3
 

## Boundary conditions

- Inlet (left): imposed parabolic axial velocity, no-slip tangential.
- Outlet (right): fixed pressure dirichlet B.C
- Embedded wall: no-slip (both normal and tangential velocity components).


## Compilation & execution

```bash
qcc -g -Wall -O2 -fopenmp -disable-dimensions DropChannel.c -o run -L$HOME/basilisk/src/gl -lglutils -lfb_tiny -lm

export OMP_NUM_THREADS=12   

>adjust to the machine being used

./run 

```


To validate against a different case from Nath et al. (2017), adjust RE, CA,
WI, MU_RATIO, BETA and DROPLET_RADIUS as needed and recompile/rerun.

## Reference

[1] B. Nath, G. Biswas, A. Dalal, and K. C. Sahu, "Migration of a droplet in a
cylindrical tube in the creeping flow regime," *Physical Review E*, 95,
033110, 2017. doi:10.1103/PhysRevE.95.033110.