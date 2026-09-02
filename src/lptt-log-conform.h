// Linear PTT Model-- adapted from the original log-conform.h 

(const) scalar lambda[] = 1.; // Polymer relaxation time
(const) scalar mup[] = 1.; // Polymer viscosity
double epsPTT = 0.05; // Linear PTT extensibility parameter

// f_s (stress) and f_r (relaxation)
void (* f_s) (double, double *, double *) = NULL;
void (* f_r) (double, double *, double *) = NULL;

// Linear PTT stress function: f_s(trA) = 1, f_r(trA) = 1 + epsPTT*(trA - d)
void lptt_f_s (double trA, double * nu, double * eta) {
  *nu  = 1.;
  *eta = 1.;
}

// eta(trA) = 1+epsPTT*(trA-d)
void lptt_f_r (double trA, double * nu, double * eta) {
  double d = dimension;
#if AXI
  d = 3.; 
#endif
  *nu  = 1.;
  *eta = 1. + epsPTT * (trA - d);
}

#include "bcg.h"

symmetric tensor tau_p[]; // Polymeric stress tensor
#if AXI
scalar tau_qq[]; // Circumferential stress component (axisymmetric)
#endif
(const) scalar trA[] = 0.; // Trace of the conformation tensor

event defaults (i = 0) {
  if (is_constant (a.x))
    a = new face vector;

  // Set default model functions to Linear PTT
  if (!f_s) f_s = lptt_f_s;
  if (!f_r) f_r = lptt_f_r;

  if (is_constant (trA))
    trA = new scalar;

  // Initialize stress fields to zero
  foreach() {
    foreach_dimension()
      tau_p.x.x[] = 0.;
      tau_p.x.y[] = 0.;
#if AXI
    tau_qq[] = 0.;
#endif
  }

  // Set default boundary conditions for stress tensor (Neumann)
  for (scalar s in {tau_p}) {
    s.v.x.i = -1; 
    foreach_dimension() {
      if (s.boundary[left] != periodic_bc) {
        s[left]  = neumann(0);
        s[right] = neumann(0);
      }
    }
  }
#if AXI
  scalar s = tau_p.x.y;
  s[bottom] = dirichlet(0.); 
#endif  
}

typedef struct { double x, y; }   pseudo_v;
typedef struct { pseudo_v x, y; } pseudo_t;

// 2D symmetric tensor A into eigenvalues (Lambda) and eigenvectors (R)
static void diagonalization_2D (pseudo_v * Lambda, pseudo_t * R, pseudo_t * A) {
  if (sq(A->x.y) < 1e-15) {
    R->x.x = R->y.y = 1.;
    R->y.x = R->x.y = 0.;
    Lambda->x = A->x.x; 
    Lambda->y = A->y.y;
    return;
  }

  double T = A->x.x + A->y.y; // Trace
  double D = A->x.x * A->y.y - sq(A->x.y); // Determinant

  R->x.x = R->x.y = A->x.y;
  R->y.x = R->y.y = -A->x.x;
  double s = 1.;
  for (int i = 0; i < dimension; i++) {
    double * ev = (double *) Lambda;
    ev[i] = T/2. + s * sqrt(sq(T)/4. - D);
    s *= -1.;
    double * Rx = (double *) &R->x;
    double * Ry = (double *) &R->y;
    Ry[i] += ev[i];
    double mod = sqrt(sq(Rx[i]) + sq(Ry[i]));
    Rx[i] /= mod;
    Ry[i] /= mod;
  }
}

event tracer_advection (i++) {
  tensor Psi = tau_p;
#if AXI
  scalar Psiqq = tau_qq;
#endif

  //tau_p to Conformation Tensor A to Psi = log(A) 
  foreach() {
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
      if (f_s) f_s (trA[], &nu, &eta);

      double fa = (mup[] != 0. ? lambda[] / (mup[] * eta) : 0.);

      // Reconstruct conformation tensor A from stress tau_p
      pseudo_t A;
      A.x.y = fa * tau_p.x.y[] / nu;
      foreach_dimension()
        A.x.x = (fa * tau_p.x.x[] + 1.) / nu;

#if AXI
      double Aqq = (1. + fa * tau_qq[]) / nu;
      Psiqq[] = log (Aqq); 
#endif

      // Compute matrix logarithm Psi = R*log(Lambda)*R^T
      pseudo_v Lambda;
      pseudo_t R;
      diagonalization_2D (&Lambda, &R, &A);
      
      Psi.x.y[] = R.x.x * R.y.x * log(Lambda.x) + R.y.y * R.x.y * log(Lambda.y);
      foreach_dimension()
        Psi.x.x[] = sq(R.x.x) * log(Lambda.x) + sq(R.x.y) * log(Lambda.y);
      
      // Compute tensor B and rotation tensor OM
      pseudo_t B;
      double OM = 0.;
      if (fabs(Lambda.x - Lambda.y) <= 1e-20) {
        B.x.y = (u.y[1,0] - u.y[-1,0] + u.x[0,1] - u.x[0,-1]) / (4. * Delta); 
        foreach_dimension() 
          B.x.x = (u.x[1,0] - u.x[-1,0]) / (2. * Delta);
      }
      else {
        pseudo_t M;
        foreach_dimension() {
          M.x.x = (sq(R.x.x)*(u.x[1] - u.x[-1]) +
                   sq(R.y.x)*(u.y[0,1] - u.y[0,-1]) +
                   R.x.x*R.y.x*(u.x[0,1] - u.x[0,-1] + u.y[1] - u.y[-1])) / (2. * Delta);
          M.x.y = (R.x.x*R.x.y*(u.x[1] - u.x[-1]) + 
                   R.x.y*R.y.x*(u.y[1] - u.y[-1]) +
                   R.x.x*R.y.y*(u.x[0,1] - u.x[0,-1]) +
                   R.y.x*R.y.y*(u.y[0,1] - u.y[0,-1])) / (2. * Delta);
        }
        double omega = (Lambda.y * M.x.y + Lambda.x * M.y.x) / (Lambda.y - Lambda.x);
        OM = (R.x.x * R.y.y - R.x.y * R.y.x) * omega;
        
        B.x.y = M.x.x * R.x.x * R.y.x + M.y.y * R.y.y * R.x.y;
        foreach_dimension()
          B.x.x = M.x.x * sq(R.x.x) + M.y.y * sq(R.x.y);  
      }

      // Advance Psi in time 
      double s = -Psi.x.y[];
      Psi.x.y[] += dt * (2. * B.x.y + OM * (Psi.y.y[] - Psi.x.x[]));
      foreach_dimension() {
        s *= -1.;
        Psi.x.x[] += dt * 2. * (B.x.x + s * OM);
      }

#if AXI
      Psiqq[] += dt * 2. * u.y[] / y;
#endif
    }
  }
  
  // Advect log-conformation field "Psi"
#if AXI
  advection ({Psi.x.x, Psi.x.y, Psi.y.y, Psiqq}, uf, dt);
#else
  advection ({Psi.x.x, Psi.x.y, Psi.y.y}, uf, dt);
#endif
  
  // Psi to A to tau_p with analytical relaxation 
  foreach() {
    if (lambda[] == 0.) {
      // Pure Newtonian limit 
      foreach_dimension()
        tau_p.x.x[] = mup[] * (u.x[1,0] - u.x[-1,0]) / Delta; 
      tau_p.x.y[]   = mup[] * (u.y[1,0] - u.y[-1,0] + u.x[0,1] - u.x[0,-1]) / (2. * Delta); 
#if AXI
      tau_qq[] = 2. * mup[] * u.y[] / y;
#endif
    }
    else { 
      // Reconstruct matrix A = exp(Psi)
      pseudo_t A = {{Psi.x.x[], Psi.x.y[]}, {Psi.y.x[], Psi.y.y[]}}, R;
      pseudo_v Lambda;
      diagonalization_2D (&Lambda, &R, &A);
      Lambda.x = exp(Lambda.x);
      Lambda.y = exp(Lambda.y);
      
      A.x.y = R.x.x * R.y.x * Lambda.x + R.y.y * R.x.y * Lambda.y;
      foreach_dimension()
        A.x.x = sq(R.x.x) * Lambda.x + sq(R.x.y) * Lambda.y;
#if AXI
      double Aqq = exp(Psiqq[]);
#endif

      double eta = 1., nu = 1.;
      if (f_r) f_r (trA[], &nu, &eta);

      // Exact integration factor for relaxation
      double fa = exp(-nu * eta * dt / lambda[]);

#if AXI
      Aqq = (1. - fa) / nu + fa * exp(Psiqq[]);
      Psiqq[] = log (Aqq);
#endif

      A.x.y *= fa;
      foreach_dimension()
        A.x.x = (1. - fa) / nu + A.x.x * fa;
      
      // Update trace of conformation tensor for non-linear models
      if (f_s || f_r) {
        scalar t = trA;
        t[] = A.x.x + A.y.y;
#if AXI
        t[] += Aqq;
#endif
      }

      nu = 1.; eta = 1.;
      if (f_s) f_s (trA[], &nu, &eta);

      // Convert conformation tensor back to polymer stress tensor tau_p
      fa = mup[] / lambda[] * eta;
      tau_p.x.y[] = fa * nu * A.x.y;
#if AXI
      tau_qq[] = fa * (nu * Aqq - 1.);
#endif
      foreach_dimension()
        tau_p.x.x[] = fa * (nu * A.x.x - 1.);
    }
  }
}

// Compute viscoelastic force divergence and add it to face acceleration field 'a'
event acceleration (i++) {
  face vector av = a;
  foreach_face() {
    if (fm.x[] > 1e-20) {
      double shear = (tau_p.x.y[0,1] * cm[0,1] + tau_p.x.y[-1,1] * cm[-1,1] -
                      tau_p.x.y[0,-1] * cm[0,-1] - tau_p.x.y[-1,-1] * cm[-1,-1]) / 4.;
      av.x[] += (shear + cm[] * tau_p.x.x[] - cm[-1] * tau_p.x.x[-1]) *
                alpha.x[] / (sq(fm.x[]) * Delta);
    }
  }
#if AXI
  foreach_face(y) {
    if (y > 0.)
      av.y[] -= (tau_qq[] + tau_qq[0,-1]) * alpha.y[] / sq(y) / 2.;
  }
#endif
}