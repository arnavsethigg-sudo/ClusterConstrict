#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "view.h"
#include "tag.h"
#include "Interface_points.h"

scalar droplet[];
face vector av[];

// symmetric about y
u.n[bottom] = dirichlet(0);
u.t[bottom] = neumann(0);
u.n[left]  = neumann(0);
u.n[right] = neumann(0);
u.n[top]   = neumann(0);
av.n[bottom] = dirichlet(0);
av.t[bottom] = dirichlet(0);
av.n[top]    = dirichlet(0);
av.t[top]    = dirichlet(0);
av.n[left]   = dirichlet(0);
av.t[left]   = dirichlet(0);
av.n[right]  = dirichlet(0);
av.t[right]  = dirichlet(0);

int main()
{
  a = av;

  size (4.);
  origin (-L0/2., 0.);
  init_grid (1 << 8);
  TOLERANCE = 1e-6;
  rho1 = rho2 = 1.;
  mu1 = mu2 = 0.01;
  f.sigma = 12.;

  run();
}

event init (t = 0)
{
  fraction (f,
            max(-(sq(x + 0.6) + sq(y) - sq(0.5)),
                -(sq(x - 0.6) + sq(y) - sq(0.5))));

  foreach()
    u.x[] = -0.5 * sign(x) * f[];
}

event compute_interface_points (i++)
{
    foreach()
        droplet[] = (f[] > 1e-6);

    tag(droplet);

    build_interface_points(f, droplet);

    set_lubrication_lengthscales (L0 / (1 << grid->maxdepth));

    double cutoff = lub_damping_cutoff;
    find_neighbour_pairs(cutoff);
    compute_minimum_gaps();
    prepare_lubrication_forces();

    printf("t=%g  h_min=%g  Fmag=%g  dir.y=%g\n",
       t, droplet_pair_gaps[0].h_min,
       n_lub_forces > 0 ? lub_forces[0].Fmag : 0.,
       n_lub_forces > 0 ? lub_forces[0].dir.y : 0.);
}

event acceleration (i++)
{
  
  foreach_face() {
    av.x[] = 0.;
    av.y[] = 0.;
  }

  apply_lubrication_forces (a);
}

event movie (t += 0.02)
{
  clear();

  draw_vof("f", lw = 1);

  mirror ({0, 1}) {
    draw_vof("f", lw = 1);
  }

  save("demonstration.mp4");
}

event end (t = 20.)
{
  printf("Interface points = %d  t = %g\n", n_interface_points, t);
}