/**
# Three-phase interfacial flows: here, at least one of the phases (f1*(1-f2)) forms a precursor film

# Version 0.2
# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids
# Last Updated: Jul 24, 2024

This file helps setup simulations for flows of three fluids separated by
corresponding interfaces (i.e. immiscible fluids). It is typically used in
combination with a [Navier--Stokes solver](navier-stokes/centered.h).

The interface between the fluids is tracked with a Volume-Of-Fluid
method. The volume fraction in fluid i is $f_i=1$ and $f_i=0$ everywhere else.
The densities and dynamic viscosities for fluid 1, 2 and 3 are *rho1*,
*mu1*, *rho2*, *mu2*, and *rho3*, *mu3* respectively. */

#include "vof.h"
scalar f1[], f2[], *interfaces = {f1, f2};
(const) scalar Gp = unity; // elastic modulus
(const) scalar lambda = unity; // relaxation time

double rho1 = 1., mu1 = 0., rho2 = 1., mu2 = 0., rho3 = 1., mu3 = 0.;
double G1 = 0., G2 = 0., G3 = 0.; // elastic moduli
double lambda1 = 0., lambda2 = 0., lambda3 = 0.; // relaxation times
double TOLelastic = 1e-1; // tolerance for elastic modulus

/**
Auxilliary fields are necessary to define the (variable) specific
volume $\alpha=1/\rho$ as well as the cell-centered density. */

face vector alphav[];
scalar rhov[];
scalar Gpd[];
scalar lambdapd[];

event defaults (i = 0) {
  alpha = alphav;
  rho = rhov;
  Gp = Gpd;
  lambda = lambdapd;
  /**
  If the viscosity is non-zero, we need to allocate the face-centered
  viscosity field. */
  mu = new face vector;
}

/**
The density and viscosity are defined using arithmetic averages by
default. The user can overload these definitions to use other types of
averages (i.e. harmonic). */

#ifndef rho
#define rho(f1, f2) (clamp(f1*(1-f2), 0., 1.) * rho1 + clamp(f1*f2, 0., 1.) * rho2 + clamp((1-f1), 0., 1.) * rho3)
#endif
#ifndef mu
#define mu(f1, f2) (clamp(f1*(1-f2), 0., 1.) * mu1 + clamp(f1*f2, 0., 1.) * mu2 + clamp((1-f1), 0., 1.) * mu3)
#endif

/**
We have the option of using some "smearing" of the density/viscosity
jump. */

#ifdef FILTERED
scalar sf1[], sf2[], *smearInterfaces = {sf1, sf2};
#else
#define sf1 f1
#define sf2 f2
scalar *smearInterfaces = {sf1, sf2};
#endif

event tracer_advection (i++)
{
  if (i > 1){
  foreach(){
    if ((f2[] > 1e-2) && (f1[] < 1.-1e-2)){
      f1[] = f2[];
    }
  }
}

/**
When using smearing of the density jump, we initialise *sfi* with the
vertex-average of *fi*. */
#ifdef FILTERED
  int counter1 = 0;
  for (scalar sf in smearInterfaces){
    counter1++;
    int counter2 = 0;
    for (scalar f in interfaces){
      counter2++;
      if (counter1 == counter2){
        // fprintf(ferr, "%s %s\n", sf.name, f.name);
      #if dimension <= 2
          foreach(){
            sf[] = (4.*f[] +
        	    2.*(f[0,1] + f[0,-1] + f[1,0] + f[-1,0]) +
        	    f[-1,-1] + f[1,-1] + f[1,1] + f[-1,1])/16.;
          }
      #else // dimension == 3
          foreach(){
            sf[] = (8.*f[] +
        	    4.*(f[-1] + f[1] + f[0,1] + f[0,-1] + f[0,0,1] + f[0,0,-1]) +
        	    2.*(f[-1,1] + f[-1,0,1] + f[-1,0,-1] + f[-1,-1] +
        		f[0,1,1] + f[0,1,-1] + f[0,-1,1] + f[0,-1,-1] +
        		f[1,1] + f[1,0,1] + f[1,-1] + f[1,0,-1]) +
        	    f[1,-1,1] + f[-1,1,1] + f[-1,1,-1] + f[1,1,1] +
        	    f[1,1,-1] + f[-1,-1,-1] + f[1,-1,-1] + f[-1,-1,1])/64.;
          }
      #endif
      }
    }
  }
#endif

#if TREE
  for (scalar sf in smearInterfaces){
    sf.prolongation = refine_bilinear;
    sf.dirty = true; // boundary conditions need to be updated
  }
#endif
}


event properties (i++) {
  
  foreach_face() {
  double ff1 = (sf1[] + sf1[-1])/2.;
  double ff2 = (sf2[] + sf2[-1])/2.;
  alphav.x[] = fm.x[]/rho(ff1, ff2);
  face vector muv = mu;
  muv.x[] = fm.x[]*mu(ff1, ff2);
  }

  foreach(){
    rhov[] = cm[]*rho(sf1[], sf2[]);

    Gpd[] = 0.;
    lambdapd[] = 0.;

    if (clamp(sf1[]*(1-sf2[]), 0., 1.) > TOLelastic){
      Gpd[] += G1*clamp(sf1[]*(1-sf2[]), 0., 1.);
      lambdapd[] += lambda1*clamp(sf1[]*(1-sf2[]), 0., 1.);
    }
    if (clamp(sf1[]*sf2[], 0., 1.) > TOLelastic){
      Gpd[] += G2*clamp(sf1[]*sf2[], 0., 1.);
      lambdapd[] += lambda2*clamp(sf1[]*sf2[], 0., 1.);
    }
    if (clamp((1-sf1[]), 0., 1.) > TOLelastic){
      Gpd[] += G3*clamp((1-sf1[]), 0., 1.);
      lambdapd[] += lambda3*clamp((1-sf1[]), 0., 1.);
    }
  
  }

#if TREE
  for (scalar sf in smearInterfaces){
    sf.prolongation = fraction_refine;
    sf.dirty = true;
  }
#endif
}