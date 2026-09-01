(const) scalar lambda[] = 1.;
(const) scalar mup[] = 1.;
double epsPTT = 0.05;

// Constitutive model functions.
void (* f_s) (double, double *, double *) = NULL;
void (* f_r) (double, double *, double *) = NULL;

// LPTT stress function.
void lptt_f_s (double trA, double * nu, double * eta) {
  *nu = 1.;
  *eta = 1.;
}

// LPTT relaxation function.
void lptt_f_r (double trA, double * nu, double * eta) {
  double d = dimension;
#if AXI
  d = 3.; 
#endif
  *nu = 1.;
  *eta = 1. + epsPTT * (trA - d);
}

// BCG advection.
#include "bcg.h"

// Polymer stress tensor.
symmetric tensor tau_p[];
#if AXI
scalar tau_qq[];
#endif
(const) scalar trA[] = 0.;

// Set initial values.
event defaults (i = 0) {
  if (is_constant (a.x))
    a = new face vector;
    
  if (!f_s) f_s = lptt_f_s;
  if (!f_r) f_r = lptt_f_r;

  if (is_constant(trA))
    trA = new scalar;

  foreach() {
    foreach_dimension()
      tau_p.x.x[] = 0.;
    tau_p.x.y[] = 0.;
#if AXI
    tau_qq[] = 0;
#endif
  }

  for (scalar s in {tau_p}) {
    s.v.x.i = -1; // just a scalar, not the component of a vector
    foreach_dimension()
      if (s.boundary[left] != periodic_bc) {
	s[left] = neumann(0);
	s[right] = neumann(0);
      }
  }
#if AXI
  scalar s = tau_p.x.y;
  s[bottom] = dirichlet (0.);  
#endif  
}

typedef struct { double x, y;}   pseudo_v;
typedef struct { pseudo_v x, y;} pseudo_t;

// Diagonalise the tensor.
static void diagonalization_2D (pseudo_v * Lambda, pseudo_t * R, pseudo_t * A)
{
  if (sq(A->x.y) < 1e-15) {
    R->x.x = R->y.y = 1.;
    R->y.x = R->x.y = 0.;
    Lambda->x = A->x.x; Lambda->y = A->y.y;
    return;
  }

  double T = A->x.x + A->y.y; 
  double D = A->x.x*A->y.y - sq(A->x.y); 

  R->x.x = R->x.y = A->x.y;
  R->y.x = R->y.y = -A->x.x;
  double s = 1.;
  for (int i = 0; i < dimension; i++) {
    double * ev = (double *) Lambda;
    ev[i] = T/2 + s*sqrt(sq(T)/4. - D);
    s *= -1;
    double * Rx = (double *) &R->x;
    double * Ry = (double *) &R->y;
    Ry[i] += ev[i];
    double mod = sqrt(sq(Rx[i]) + sq(Ry[i]));
    Rx[i] /= mod;
    Ry[i] /= mod;
  }
}

// Advect and update the stress.
event tracer_advection (i++)
{

// Use tau_p for log-conformation.
  tensor Psi = tau_p;
#if AXI
  scalar Psiqq = tau_qq;
#endif

  foreach() {
    if (cs[] <= 0.)
      continue;
    if (lambda[] == 0.) {
      foreach_dimension()
	Psi.x.x[] = 0.;
      Psi.x.y[] = 0.;
#if AXI
      Psiqq[] = 0.;
#endif
    }
    else { 

      double eta = 1., nu = 1.;
      if (f_s)
	f_s (trA[], &nu, &eta);

// Convert stress to A.
      double fa = (mup[] != 0 ? lambda[]/(mup[]*eta) : 0.);

      pseudo_t A;
      A.x.y = fa*tau_p.x.y[]/nu;
      foreach_dimension()
	A.x.x = (fa*tau_p.x.x[] + 1.)/nu;

#if AXI
      double Aqq = (1. + fa*tau_qq[])/nu;

      if (Aqq <= 0. || !isfinite(Aqq)) {
        fprintf(stderr,
                "\nBAD AQQ BEFORE LOG\n"
                "t=%g x=%g y=%g\n"
                "Aqq=%g tauqq=%g\n"
                "fa=%g nu=%g\n"
                "cs=%g lambda=%g mup=%g\n",
                t,x,y,
                Aqq,tau_qq[],
                fa,nu,
                cs[],lambda[],mup[]);
        fflush(stderr);
      }

      Psiqq[] = log(Aqq);
#endif

// Eigenvalues and eigenvectors.
      pseudo_v Lambda;
      pseudo_t R;
      diagonalization_2D (&Lambda, &R, &A);
      
double lx = log(max(Lambda.x,1e-30));
double ly = log(max(Lambda.y,1e-30));

Psi.x.y[] = R.x.x*R.y.x*lx + R.y.y*R.x.y*ly;

foreach_dimension()
  Psi.x.x[] = sq(R.x.x)*lx + sq(R.x.y)*ly;
// Velocity-gradient terms.
      pseudo_t B;
      double OM = 0.;
      if (fabs(Lambda.x - Lambda.y) <= 1e-20) {
	B.x.y = (u.y[1,0] - u.y[-1,0] +
		 u.x[0,1] - u.x[0,-1])/(4.*Delta); 
	foreach_dimension() 
	  B.x.x = (u.x[1,0] - u.x[-1,0])/(2.*Delta);
      }
      else {
	pseudo_t M;
	foreach_dimension() {
	  M.x.x = (sq(R.x.x)*(u.x[1] - u.x[-1]) +
		   sq(R.y.x)*(u.y[0,1] - u.y[0,-1]) +
		   R.x.x*R.y.x*(u.x[0,1] - u.x[0,-1] + 
				u.y[1] - u.y[-1]))/(2.*Delta);
	  M.x.y = (R.x.x*R.x.y*(u.x[1] - u.x[-1]) + 
		   R.x.y*R.y.x*(u.y[1] - u.y[-1]) +
		   R.x.x*R.y.y*(u.x[0,1] - u.x[0,-1]) +
		   R.y.x*R.y.y*(u.y[0,1] - u.y[0,-1]))/(2.*Delta);
	}
	double omega = (Lambda.y*M.x.y + Lambda.x*M.y.x)/(Lambda.y - Lambda.x);
	OM = (R.x.x*R.y.y - R.x.y*R.y.x)*omega;
	
	B.x.y = M.x.x*R.x.x*R.y.x + M.y.y*R.y.y*R.x.y;
	foreach_dimension()
	  B.x.x = M.x.x*sq(R.x.x)+M.y.y*sq(R.x.y);	
      }

      double s = - Psi.x.y[];
      Psi.x.y[] += dt*(2.*B.x.y + OM*(Psi.y.y[] - Psi.x.x[]));
      foreach_dimension() {
	s *= -1;
	Psi.x.x[] += dt*2.*(B.x.x + s*OM);
      }

#if AXI
      Psiqq[] += dt*2.*u.y[]/y;
#endif
    }
  }

#if AXI
// Advect the log-conformation.
  advection ({Psi.x.x, Psi.x.y, Psi.y.y, Psiqq}, uf, dt);
#else
  advection ({Psi.x.x, Psi.x.y, Psi.y.y}, uf, dt);
#endif

  
  foreach() {
    if (cs[] <= 0.)
      continue;

    if (lambda[] == 0.) {

      foreach_dimension()
	tau_p.x.x[] = mup[]*(u.x[1,0] - u.x[-1,0])/Delta; 
      tau_p.x.y[] = mup[]*(u.y[1,0] - u.y[-1,0] +
			   u.x[0,1] - u.x[0,-1])/(2.*Delta); 
#if AXI
      tau_qq[] = 2.*mup[]*u.y[]/y;
#endif
    }
    else { 
      
// Convert log(A) back to A.
      pseudo_t A = {{Psi.x.x[], Psi.x.y[]}, {Psi.y.x[], Psi.y.y[]}}, R;
      pseudo_v Lambda;
      diagonalization_2D (&Lambda, &R, &A);
      Lambda.x = exp(Lambda.x), Lambda.y = exp(Lambda.y);
      
      A.x.y = R.x.x*R.y.x*Lambda.x + R.y.y*R.x.y*Lambda.y;
      foreach_dimension()
	A.x.x = sq(R.x.x)*Lambda.x + sq(R.x.y)*Lambda.y;
#if AXI
      double Aqq = exp(Psiqq[]);
#endif

      double eta = 1., nu = 1.;
      if (f_r) {
	f_r (trA[], &nu, &eta);
      }

// Apply PTT relaxation.
      double fa = exp(-nu*eta*dt/lambda[]);

#if AXI
      Aqq = (1. - fa)/nu + fa*exp(Psiqq[]);
      Psiqq[] = log (Aqq);
#endif

      A.x.y *= fa;
      foreach_dimension()
	A.x.x = (1. - fa)/nu + A.x.x*fa;
      
// Update tr(A).
      if (f_s || f_r) {
	scalar t = trA;
	t[] = A.x.x + A.y.y;
#if AXI
	t[] += Aqq;
#endif
      }

      nu = 1; eta = 1.;
      if (f_s)
	f_s (trA[], &nu, &eta);

      fa = mup[]/lambda[]*eta;
      
// Convert A back to stress.
      tau_p.x.y[] = fa*nu*A.x.y;
#if AXI
      tau_qq[] = fa*(nu*Aqq - 1.);
#endif
      foreach_dimension()
	tau_p.x.x[] = fa*(nu*A.x.x - 1.);
    }
  }
 

}

// Add polymer stress to momentum.
event acceleration (i++)
{
  face vector av = a;
  foreach_face()
    if (fm.x[] > 1e-20) {
      double shear = (tau_p.x.y[0,1]*cm[0,1] + tau_p.x.y[-1,1]*cm[-1,1] -
		      tau_p.x.y[0,-1]*cm[0,-1] - tau_p.x.y[-1,-1]*cm[-1,-1])/4.;
      av.x[] += (shear + cm[]*tau_p.x.x[] - cm[-1]*tau_p.x.x[-1])*
	alpha.x[]/(sq(fm.x[])*Delta);
    }
#if AXI
  foreach_face(y)
    if (y > 0.)
      av.y[] -= (tau_qq[] + tau_qq[0,-1])*alpha.y[]/sq(y)/2.;
#endif
}