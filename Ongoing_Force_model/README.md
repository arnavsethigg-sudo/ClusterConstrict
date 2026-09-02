# Two-Droplet Lubrication/Repulsion Force Model

A Basilisk C solver with a
pairwise short-range force model, combining a hydrodynamic-film lubrication
repulsion with a Morse-potential term.
## Problem setup

- 2-D, uniform grid (grid/multigrid.h), domain size 4, resolution 1<<8,
  symmetric about y = 0 (free-slip bottom boundary, no normal flow through the
  other three sides).
- Two droplets of radius 0.5 are initialised at x = -0.6 and x = 0.6, each
  given an initial velocity toward the other (u.x = -0.5 * sign(x) * f),
  so they approach and interact head-on.
- Equal density, equalviscosity (rho1 = rho2 = 1, mu1 = mu2 = 0.01),
  surface tension f.sigma = 12.

## Brief Model explanation

1. Tag connected droplet regions (tag()) so each droplet has a distinct
   integer ID.
2. Extract interface points: for every interface (mixed) cell, reconstruct
   the PLIC plane and record its centroid, normal, tag, grid level and
   local Delta (build_interface_points in Interface_points.h).
3. Find candidate neighbour pairs: search over all interface points,
   keeping pairs belonging to different droplets within a cutoff distance
   (find_neighbour_pairs).
4. Reduce to one minimum gap per droplet pair (compute_minimum_gaps),
   tracking which interface points satisfies that minimum.
5. Convert gaps to forces (prepare_lubrication_forces): each pair's force
   magnitude is the sum of two closed-form terms, both functions of the gap h —
   - lubrication_force_magnitude — a short-range 1/h³ repulsion.
     between near-touching interfaces.
   - morse_force_magnitude — a Morse-potential force
     with an equilibrium separation morse_r0 and interaction cutoff
     morse_cutoff.
6. Apply the resulting forces (apply_lubrication_forces, called from the
   acceleration event): each pair's force is spread onto nearby velocity
   faces with a normalized Gaussian gaussian spread (inject_face_force) rather than
   applied as a point force, and mirrored across y = 0 to respect the
   symmetric half-domain setup.

All lubrication/Morse length scales (lub_h0, lub_hcutoff, morse_r0,
morse_cutoff) are set from a common grid length scale via
set_lubrication_lengthscales, called once per step with the finest grid
spacing L0 / (1 << grid->maxdepth).

## Compilation & execution


```bash
qcc -g -Wall -O2 -fopenmp -disable-dimensions DropChannel.c     -o run     -L$HOME/basilisk/src/gl     -lglutils -lfb_tiny -lm

./run

```



