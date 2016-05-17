#include "mp.h"
#include "math.h"
extern real ScalingFactor;

/* Surface density */
real Sigma(r)
     real r;
{
  real Sig;
  real cavity=1.0;
  if (r < CAVITYRADIUS) cavity = 1.0/CAVITYRATIO; /* This is *not* a steady state */
                                /* profile, if a cavity is defined. It first needs */
                                /* to relax towards steady state, on a viscous time scale */
  if ( fabs(RATIOTRAN-1.0) > 0.001 && r < RTRAN+DRTRAN) {
        if (r < RTRAN-DRTRAN) cavity = pow(1.0/RATIOTRAN,0.8);
        if (r >= RTRAN-DRTRAN) cavity = pow((2.0/(RATIOTRAN*1.0-1.0))/(sin(-1.0*(r-RTRAN)*M_PI/(2.0*DRTRAN))+(1.0*RATIOTRAN+1.0)/(1.0*RATIOTRAN-1.0)),0.8);
        }
  Sig=cavity*ScalingFactor*SIGMA0*pow(r,-SIGMASLOPE);
	
  return Sig;

}

void FillSigma() {
  int i;
  for (i = 0; i < NRAD; i++) {
    SigmaMed[i] = Sigma(Rmed[i]);
    SigmaInf[i] = Sigma(Rinf[i]);
  }
}

void RefillSigma (Surfdens)
     PolarGrid *Surfdens;
{
  int i, j, nr, ns, l;
  real *field;
  real moy;
  nr = Surfdens->Nrad;
  ns = Surfdens->Nsec;
  field = Surfdens->Field;
  for (i = 0; i < nr; i++) {
    moy = 0.0;
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      moy += field[l];
    }
    moy /= (real)ns;
    SigmaMed[i] = moy;
SigmaMed1[i]=moy;
  }
  SigmaInf[0] = SigmaMed[0];
  for (i = 1; i < nr; i++) {
    SigmaInf[i] = (SigmaMed[i-1]*(Rmed[i]-Rinf[i])+\
		   SigmaMed[i]*(Rinf[i]-Rmed[i-1]))/\
      (Rmed[i]-Rmed[i-1]);
  }
}

/* Thermal energy */
real Energy(r)
     real r;
{
  real energy0;
  FILE* fp=fopen("energy_1AU.txt","w+");
  if (ADIABATICINDEX == 1.0) {
    fprintf (stderr, "The adiabatic index must differ from unity to initialize the gas internal energy. I must exit.\n");
    prs_exit (1);
  }
  else
    //energy0 = R/MU/(ADIABATICINDEX-1.0)*SIGMA0*pow(ASPECTRATIO,2.0)*pow(r,-SIGMASLOPE-1.0+2.0*FLARINGINDEX);
 // energy0 = R/MU/(ADIABATICINDEX-1.0)*Sigma(r)*pow(ASPECTRATIO,2.0)*pow(r,-1.0+2.0*FLARINGINDEX);
  energy0 = 1.0/(ADIABATICINDEX-1.0)*Sigma(r)*pow(ASPECTRATIO,2.0)*pow(r,-1.0+2.0*FLARINGINDEX);

  //if(fabs(r-1.0)<0.1 || 1) printf("%f\t%e\n",r,energy0);
  if(energy0<0.0) printf("r=%f\t E=%g\n",r,energy0);
  if(fabs(r-1.0)<0.1) fprintf(fp,"%f\n",energy0);
  fclose(fp);
  return energy0;
}

void FillEnergy() {
  int i;
  for (i = 0; i < NRAD; i++)
    {
	EnergyMed[i] = Energy(Rmed[i]);
//	if (i<10) EnergyMed[i]=EnergyMed[i]/2.0;
}
}


void RefillEnergy (energy)
     PolarGrid *energy;
{
  int i, j, nr, ns, l;
FILE *fp;
char name[256];
sprintf (name, "%sEnergyMed%d.dat", OUTPUTDIR,CPU_Rank);
fp=fopen(name,"w");
 
real *field;
  real moy;
  nr = energy->Nrad;
  ns = energy->Nsec;
  field = energy->Field;
  for (i = 0; i < nr; i++) {
    moy = 0.0;
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      moy += field[l];
    }
    moy /= (real)ns;
    EnergyMed[i] = moy;///(ADIABATICINDEX-1.0)*SigmaMed[i];
	EnergyMed1[i]=SigmaMed[i]*moy/(ADIABATICINDEX-1.0);
	fprintf(fp,"%2.12g\t%2.12g\n",EnergyMed[i]*(ADIABATICINDEX-1.0)/SigmaMed[i]*TEUNIT,EnergyMed1[i]);
  }
fclose(fp);
}
/*void RefillVtheta (vtheta)
     PolarGrid *vtheta;
{
  int i, j, nr, ns, l;
  real *field;
  real moy;
  nr = vtheta->Nrad;
  ns = vtheta->Nsec;
  field = vtheta->Field;
  for (i = 0; i < nr; i++) {
    moy = 0.0;
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      moy += field[l];
    }
    moy /= (real)ns;
    VthetaMed[i] = moy;
  }
}

void RefillVrad (vrad)
     PolarGrid *vrad;
{
  int i, j, nr, ns, l;
  real *field;
  real moy;
  nr = vrad->Nrad;
  ns = vrad->Nsec;
  field = vrad->Field;
  for (i = 0; i < nr; i++) {
    moy = 0.0;
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      moy += field[l];
    }
    moy /= (real)ns;
    VradMed[i] = moy;
  }
}
*/
/* Cooling time */
real CoolingTime(r)
     real r;
{
  real ct0;
  ct0 = COOLINGTIME0*pow(r,2.0+2.0*FLARINGINDEX);
  return ct0;
}

void FillCoolingTime() {
  int i;
  for (i = 0; i < NRAD; i++)
    CoolingTimeMed[i] = CoolingTime(Rmed[i]);
}

/* Heating source term */
real Qplusinit(r)
     real r;
{
  real qp0, viscosity;
  viscosity = FViscosity(r);
  qp0 = 2.25*viscosity*SIGMA0*pow(r,-SIGMASLOPE-3.0);
  return qp0;
}

void FillQplus() {
  int i;
  for (i = 0; i < NRAD; i++)
    QplusMed[i] = Qplusinit(Rmed[i]);
}
