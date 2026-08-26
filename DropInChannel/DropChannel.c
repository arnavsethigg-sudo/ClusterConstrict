#include "embed.h"
#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "lptt-log-conform-embed.h"
#include "view.h"

// Non-dimensional control parameters 
#define RE        1.0   // Continuous-phase Reynolds number (rho2 * U * R / mu2) [u=1]
#define CA        1.  // Capillary number (mu2 * U / sigma)
#define WI        0.   // Weissenberg number (lambda * U / R)
#define MU_RATIO  1.0   // Viscosity ratio (solvent to polymeric)
#define BETA      0.1   // Solvent viscosity fraction
#define RHO_RATIO 1.0   // Density ratio (rho1 / rho2)

#define RIN 1.0
/* Define any curved geometry here
#define RTHROAT 0.5
#define XSTART -3.
#define XEND 5.
*/
#define DROPLET_RADIUS 0.8
#define X0 -8.


scalar lambdaf[], mupf[];
scalar cmv[];
face vector fmv[];

double LAM, MUPP;

// Straight channel wall, for validation against Nath et al. (2017)
static inline double wall_radius (double x)
{
  return RIN;
}

// Set up the axisymmetric domain, fluid properties, and solver parameters.
int main()
{
  cm = cmv;
  fm = fmv;

  size (20.);
  origin (-10.,0.);

  init_grid (1<<9);


  // currently, the density of the droplet and continuous phase are equal.(neutral buoyancy)
  rho2 = 1.0;
  rho1 = RHO_RATIO * rho2;

  // Continuous-phase viscosity 
  mu2 = (rho2 * 1.0 * RIN) / RE;

  // Viscosities derived from MU_RATIO and solvent fraction BETA
  double mu1_total = MU_RATIO * mu2;
  mu1  = BETA * mu1_total;
  MUPP = (1.0 - BETA) * mu1_total;

  // Surface tension from Capillary number
  f.sigma = (mu2 * 1.0) / CA;

  // Relaxation time from Weissenberg number
  LAM = (WI * RIN) / 1.0;

  // Linear PTT extensibility parameter.
  epsPTT = 0.25;

  CFL = 0.4;
  DT = 1e-3;
  TOLERANCE = 1e-4;
  NITERMAX = 100;

  run();
}

// Inlet prescribed parabolic axial velocity profile.
u.n[left] = dirichlet (2.*(1. - sq(y/RIN)));
u.t[left] = dirichlet (0.);
// Outlet fixed pressure with zero normal velocity gradient.
p[right]   = dirichlet (0.);
u.n[right] = neumann (0.);
u.t[right] = neumann (0.);
// No-slip condition on the embedded wall.
u.n[embed] = dirichlet (0.);
u.t[embed] = dirichlet (0.);


event init (t = 0)
  {

vertex scalar phi[];
    foreach_vertex()
      phi[] = wall_radius(x) - y;

    fractions (phi, cs, fs);
    fractions_cleanup (cs, fs);
    cm_update (cm, cs, fs);
    fm_update (fm, cs, fs);


    cm.refine = cm.prolongation = refine_cm_axi;
    cs.refine = cs.prolongation = fraction_refine;
    fm.x.refine = fm.x.prolongation = refine_face_x_axi;
    fm.y.refine = fm.y.prolongation = refine_face_y_axi;
    metric_embed_factor = axi_factor;
    boundary ({cs,fs});
    restriction ({cs,fs,cm,fm});

    lambda = lambdaf;
    mup    = mupf;

// Initialise the droplet centred at X0.
    fraction (f, sq(DROPLET_RADIUS) - (sq(x - X0) + sq(y)));

  foreach() {
// Ensuring the viscoelastic properties are only assigned to the droplet.
    lambdaf[] = LAM*f[];
    mupf[] = MUPP*f[];
    trA[] = 3.; // Tr(A) Initialisation
  }

    boundary ({f,lambdaf,mupf});
  }


  event properties (i++)
  {
    foreach() {
      lambdaf[] = LAM*f[];
      mupf[]    = MUPP*f[];
  }
  boundary ({lambdaf,mupf});
  }


  event logfile (t += 0.001)
  {
     fprintf(stderr,"t=%g \n",t);
  }

// Adapt the grid around the interface, wall, and velocity field.
event adapt (i++)
{

   adapt_wavelet ({cs, f, u.x, u.y}, (double[]){1e-2, 1e-3, 1e-3, 1e-3}, 10, 6);

   vertex scalar phi[];
   foreach_vertex()
     phi[] = wall_radius(x) - y;
   fractions (phi, cs, fs);
   fractions_cleanup (cs, fs);
   cm_update (cm, cs, fs);
   fm_update (fm, cs, fs);
   boundary ({cs,fs});
   restriction ({cs,fs,cm,fm});
}

  event movie (t += 0.01)
  {
    clear();

    view (width = 1800, height = 500, fov = 10);
    draw_vof ("cs", lw = 2);
    draw_vof ("f", lw = 2);
    mirror ({0,1}) {
    draw_vof ("cs", lw = 2);
    draw_vof ("f", filled = 1);
    }

    save ("Nath_validation.mp4");
  }

  event end (t = 10.);