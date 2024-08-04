/** Title: Bursting bubble when the coating has a tendency to dewet
# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids
# Last Updated: Jul 27, 2024
*/

// 1 is bulk (full wetting), 2 is coating (water+PEO, this will dewet) and 3 is air
#include "axi.h"
#include "navier-stokes/centered.h"
#define FILTERED
#include "three-phase-viscoelastic.h"
#include "log-conform-viscoelastic.h"
#include "tension.h"
#include "distance.h"

#define MINlevel 3                                              // maximum level

#define tsnap (5e-3)

// Error tolerances
#define fErr (1e-3)                                 // error tolerance in VOF
#define KErr (1e-4)                                 // error tolerance in KAPPA
#define VelErr (1e-2)                            // error tolerances in velocity

// boundary conditions

double Ohbulk, muR_cb, muR_ab, rho_cb, rho_ab, Ec, De, tmax, Ldomain;
int MAXlevel;

int main(int argc, char const *argv[]) {
  
  // bulk is hexadecane for cases I and II
  Ohbulk = 0.017; // this is for case I
  
  muR_cb = 1e0/3.3; // this is only the (solvent) viscosity coated to bulk viscosity ratio
  rho_cb = 1e0/0.781;
  
  muR_ab = 5e-3;
  rho_ab = 1.3e-3;

  Ec = 0.20;
  De = 0.20;

  tmax = 1e0;
  Ldomain = 8e0;
  MAXlevel = 11;

  L0=Ldomain;
  X0=-2e0; Y0=0.;
  init_grid (1 << (4));

  char comm[80];
  sprintf (comm, "mkdir -p intermediate");
  system(comm);

  rho1 = 1e0; mu1 = Ohbulk; G1 = 0.; lambda1 = 0.;
  rho2 = rho_cb; mu2 = muR_cb*Ohbulk; G2 = Ec; lambda2 = De;
  rho3 = rho_ab; mu3 = muR_ab*Ohbulk; G3 = 0.; lambda3 = 0.;

  f1.sigma = 1e0; // surface tension bulk-air!!!: 0.02747 N/m
  f2.sigma = 26.8/27.47; // bulk-coating interfacial tension: 0.0268 N/m

  fprintf(ferr, "Level %d tmax %g. Oh %3.2e, muR_cb %3.2e, rho_cb %3.2e, Ec %3.2e, De %3.2e\n", MAXlevel, tmax, Ohbulk, muR_cb, rho_cb, Ec, De);
  run();

}

event init(t = 0){
  if(!restore (file = "dump")){
    char filename1[60], filename2[60];
    /**
    Initialization for f1 and f2
    */
    sprintf(filename1, "initialconditions/CaseI_VeryThinLayer_PEO0.2wt/f1.dat");
    sprintf(filename2, "initialconditions/CaseI_VeryThinLayer_PEO0.2wt/f2.dat");
    // sprintf(filename1, "initialconditions/CaseII_thickLayer_PIB_PEO0.2wt/f1.dat");
    // sprintf(filename2, "initialconditions/CaseII_thickLayer_PIB_PEO0.2wt/f2.dat");

    FILE * fp1 = fopen(filename1,"rb");
    if (fp1 == NULL){
      fprintf(ferr, "There is no file named %s\n", filename1);
      return 1;
    }
    FILE * fp2 = fopen(filename2,"rb");
    if (fp2 == NULL){
      fprintf(ferr, "There is no file named %s\n", filename2);
      return 1;
    }

    coord* InitialShape1;
    coord* InitialShape2;
    InitialShape1 = input_xy(fp1);
    fclose (fp1);
    InitialShape2 = input_xy(fp2);
    fclose (fp2);
    scalar d1[], d2[];

    distance (d1, InitialShape1);
    distance (d2, InitialShape2);

    while (adapt_wavelet ((scalar *){f1, f2, d1, d2}, (double[]){1e-8, 1e-8, 1e-8, 1e-8}, MAXlevel).nf);
    /**
    The distance function is defined at the center of each cell, we have
    to calculate the value of this function at each vertex. */
    vertex scalar phi1[], phi2[];
    foreach_vertex(){
      phi1[] = -(d1[] + d1[-1] + d1[0,-1] + d1[-1,-1])/4.;
      phi2[] = (d2[] + d2[-1] + d2[0,-1] + d2[-1,-1])/4.;
    }
    /**
    We can now initialize the volume fractions in the domain. */
    fractions (phi1, f1);
    fractions (phi2, f2);
  }
  // dump (file = "dump");
  // return 1;
}

scalar KAPPA1[], KAPPA2[];
event adapt(i++) {
  curvature(f1, KAPPA1);
  curvature(f2, KAPPA2);

  adapt_wavelet ((scalar *){f1, f2, u.x, u.y, KAPPA1, KAPPA2},
    (double[]){fErr, fErr, VelErr, VelErr, KErr, KErr},
    MAXlevel, MINlevel);
}

// Outputs
event writingFiles (t = 0; t += tsnap; t <= tmax + tsnap) {
  dump (file = "dump");
  char nameOut[80];
  sprintf (nameOut, "intermediate/snapshot-%5.4f", t);
  dump (file = nameOut);
}

event logWriting (i++) {
  double ke = 0.;
  foreach (reduction(+:ke)){
    ke += sq(Delta)*(sq(u.x[]) + sq(u.y[]))*rho(f1[],f2[]);
  }
  static FILE * fp;
  if (i == 0) {
    fprintf (ferr, "i dt t ke\n");
    fp = fopen ("log", "w");
    fprintf(fp, "Level %d tmax %g. Oh %3.2e, muR_cb %3.2e, rho_cb %3.2e, Ec %3.2e, De %3.2e\n", MAXlevel, tmax, Ohbulk, muR_cb, rho_cb, Ec, De);
    fprintf (fp, "i dt t ke\n");
    fprintf (fp, "%d %g %g %g\n", i, dt, t, ke);
    fclose(fp);
  } else {
    fp = fopen ("log", "a");
    fprintf (fp, "%d %g %g %g\n", i, dt, t, ke);
    fclose(fp);
  }
  fprintf (ferr, "%d %g %g %g\n", i, dt, t, ke);
}