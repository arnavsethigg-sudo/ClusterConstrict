#include "axi.h"
#include "navier-stokes/centered.h"
#include "lptt-log-conform.h"


#define WI   1.0 // Weissenberg number (lambda*U_N / R)
#define BETA 0.01 // Solvent viscosity fraction (mu_s/mu0)
#define U_N  (1.0/8.0) // Newtonian mean velocity reference (Pipe Flow): U_N =(a_x*R^2)/(8*mu0)=1/8

scalar lambdaf[], mupf[];
face vector muv[];

double LAMBDA, MUPP, MUS;

int main()
{
  periodic (right);

  size (1.);
  init_grid (1 << 7);

  u.t[top] = dirichlet (0.);

  double mu0 = 1.0; // Total viscosity (mu_s + mu_p = 1)

  // Solvent and polymeric viscosity components
  MUS  = BETA * mu0;
  MUPP = (1.0 - BETA) * mu0;



  // Relaxation time (derived from Weissenberg number)
  LAMBDA = (WI * 1.0) / U_N;

  DT = 1e-3;
  epsPTT = 0.25;

  run();
}

event properties (i++)
{
  foreach_face(x)
    muv.x[] = MUS*fm.x[];

  foreach_face(y)
    muv.y[] = MUS*fm.y[];

  mu = muv;
}

event init (t = 0)
{
  lambda = lambdaf;
  mup = mupf;

  foreach() {
    lambdaf[] = LAMBDA;
    mupf[] = MUPP;
    trA[] = 3.; // Tr(A) = 3 initialisation
  }

  boundary ({lambda, mup});
}

event acceleration (i++)
{
  foreach_face(x)
    a.x[] += 1.0;
}

event diagnostics (t += 0.1)
{
  fprintf (stderr, "t=%g \n", t);
}

event profile (t = 70.)
{
  
  FILE * fp;
  char name[80];
  sprintf (name, "profile_beta_%g.dat", BETA);

  fp = fopen (name, "w");

  foreach()
    fprintf (fp, "%g %g\n", y, u.x[] / U_N);

  fclose (fp);
}

event end (t = 70.); // Change appropriately to reach steady - state simulation.