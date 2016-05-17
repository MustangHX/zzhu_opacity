#include "mp.h"

#define CFLSECURITY 0.5		/* Maximum fraction of zone size */
				/* swept in one timestep */

#define CVNR 1.41       	/* Shocks are spread over CVNR zones:       */
                                /* von Neumann-Richtmyer viscosity constant */
				/* Beware of misprint in Stone and Norman's */
				/* paper : use C2^2 instead of C2           */

static PolarGrid *TemperInt;
static PolarGrid *VradNew, *VradInt;
static PolarGrid *VthetaNew, *VthetaInt;
static PolarGrid *EnergyNew, *EnergyInt;
static real timeCRASH;  
extern boolean Corotating;
extern boolean Adiabatic;
//extern boolean SelfGravity, SGZeroMode;
extern boolean ZMPlus;
real PhysicalTime=0.0, OmegaFrame, PhysicalTimeInitial;
int FirstGasStepFLAG=1;
static int AlreadyCrashed = 0, GasTimeStepsCFL;
extern boolean FastTransport, IsDisk;
extern boolean Restart;
Pair DiskOnPrimaryAcceleration;

boolean DetectCrash (array)
     PolarGrid *array;
{
  int i, j, l, nr, ns, ic, jc, lc;
  real *ptr;  
  boolean bool = NO;
  nr = array->Nrad;
  ns = array->Nsec;
  ptr= array->Field;
#pragma omp parallel for private(j,l) shared(bool)
  for (i = 0; i < nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      if (ptr[l] < 0.0){
	bool = YES;
	ic = i;
	jc = j;
	lc = l;
	//printf("energy=%g\t%g\t%d\t%d\n",ptr[l],ptr[l],i,j);
    }
  }
  if (bool == YES){
	fprintf(stdout,"crash i j value %d %d %g %s\n",ic,jc,ptr[lc],array -> Name);
  }
  return bool;
}
}
 
void FillPolar1DArrays ()
{
  FILE *input, *output;
  int i,ii;
  real drrsep;
  float temporary;
  char InputName[256], OutputName[256];
  drrsep = (RMAX-RMIN)/(real)GLOBALNRAD;
  sprintf (InputName, "%s%s", OUTPUTDIR, "radii.dat");
  sprintf (OutputName, "%s%s", OUTPUTDIR, "used_rad.dat");
  input = fopen (InputName, "r");
  if (input == NULL) {
    mastererr ("Warning : no `radii.dat' file found. Using default.\n");
    if (LogGrid == YES) {
      for (i = 0; i <= GLOBALNRAD; i++) {
	Radii[i] = RMIN*exp((real)i/(real)GLOBALNRAD*log(RMAX/RMIN));
      }
    } else {
      for (i = 0; i <= GLOBALNRAD; i++) {
	Radii[i] = RMIN+drrsep*(real)(i);
      }
    }
  } else {
    mastererr ("Reading 'radii.dat' file.\n");
    for (i = 0; i <= GLOBALNRAD; i++) {
      fscanf (input, "%f", &temporary);
      Radii[i] = (real)temporary;
    }
  }
  for (i = 0; i < GLOBALNRAD; i++) {
    GlobalRmed[i] = 2.0/3.0*(Radii[i+1]*Radii[i+1]*Radii[i+1]-Radii[i]*Radii[i]*Radii[i]);
    GlobalRmed[i] = GlobalRmed[i] / (Radii[i+1]*Radii[i+1]-Radii[i]*Radii[i]);
  }
  for (i = 0; i < NRAD; i++) {
    ii = i+IMIN;
    Rinf[i] = Radii[ii];
    Rsup[i] = Radii[ii+1];
    Rmed[i] = 2.0/3.0*(Rsup[i]*Rsup[i]*Rsup[i]-Rinf[i]*Rinf[i]*Rinf[i]);
    Rmed[i] = Rmed[i] / (Rsup[i]*Rsup[i]-Rinf[i]*Rinf[i]);
    Surf[i] = PI*(Rsup[i]*Rsup[i]-Rinf[i]*Rinf[i])/(real)NSEC;
    InvRmed[i] = 1.0/Rmed[i];
    InvSurf[i] = 1.0/Surf[i];
    InvDiffRsup[i] = 1.0/(Rsup[i]-Rinf[i]);
    InvRinf[i] = 1.0/Rinf[i];
  }
  Rinf[NRAD]=Radii[NRAD+IMIN];
  for (i = 1; i < NRAD; i++) {
    InvDiffRmed[i] = 1.0/(Rmed[i]-Rmed[i-1]);
  }
  if (CPU_Master) {
    output = fopen (OutputName, "w");
    if (output == NULL) {
      mastererr ("Can't write %s.\nProgram stopped.\n", OutputName);
      prs_exit (1);
    }
    for (i = 0; i <= GLOBALNRAD; i++) {
      fprintf (output, "%.18g\n", Radii[i]);
    }
    fclose (output);
  }
  if (input != NULL) fclose (input);
}

void InitEuler (Vr, Vt, Rho, Energy)
     PolarGrid *Vr, *Vt, *Rho, *Energy;
{
  InitTransport ();
  InitViscosity ();
  RhoStar      = CreatePolarGrid(NRAD, NSEC, "RhoStar");
  RhoInt       = CreatePolarGrid(NRAD, NSEC, "RhoInt");
  VradNew      = CreatePolarGrid(NRAD, NSEC, "VradNew");
  VradInt      = CreatePolarGrid(NRAD, NSEC, "VradInt");
  VthetaNew    = CreatePolarGrid(NRAD, NSEC, "VthetaNew");
  VthetaInt    = CreatePolarGrid(NRAD, NSEC, "VthetaInt");
  EnergyNew    = CreatePolarGrid(NRAD, NSEC, "EnergyNew");
  EnergyInt    = CreatePolarGrid(NRAD, NSEC, "EnergyInt");
  TemperInt    = CreatePolarGrid(NRAD, NSEC, "TemperInt");
  Potential    = CreatePolarGrid(NRAD, NSEC, "Potential");
  Pressure     = CreatePolarGrid(NRAD, NSEC, "Pressure");
  SoundSpeed   = CreatePolarGrid(NRAD, NSEC, "SoundSpeed");
  Temperature  = CreatePolarGrid(NRAD, NSEC, "Temperature");
  Qplus        = CreatePolarGrid(NRAD, NSEC, "Qplus");
  Qcoola       = CreatePolarGrid(NRAD, NSEC, "Qcoola");
  Opacity       = CreatePolarGrid(NRAD, NSEC, "Opacity");

  InitComputeAccel ();
  /* Rho and Energy are already initialized: cf main.c */
  ComputeSoundSpeed (Rho, Energy);
  ComputePressureField (Rho, Energy);
  if (!(Restart)){
  ComputeTemperatureField (Rho, Energy);
  InitGasVelocities (Vr, Vt);
  }
}

real min2 (a,b)
     real a,b;
{
  if (b < a) return b;
  return a;
}

real max2 (a,b)
     real a,b;
{
  if (b > a) return b;
  return a;
}


void ActualiseGas (array, newarray)
     PolarGrid *array, *newarray;
{
  int i,j,l,ns,nr;
  real *old, *new;
  nr = array->Nrad;
  ns = array->Nsec;
  old= array->Field;
  new= newarray->Field;
#pragma omp parallel for private(j,l)
  for (i = 0; i <= nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j+ns*i;
      old[l] = new[l];
    }
  }
}


void AlgoGas (force, Rho, Vrad, Vtheta, Energy, Label, sys, TimeStep)
     Force *force;
     PolarGrid *Rho, *Vrad, *Vtheta, *Energy, *Label;
     PlanetarySystem *sys;
     int TimeStep;
{
  real dt, dtemp=0.0;
  real OmegaNew, domega;
  int gastimestepcfl;
  boolean Crashed=NO;
  FirstGasStepFLAG=1;
  gastimestepcfl = 1;
  if (Adiabatic) {
    ComputeSoundSpeed (Rho, Energy);
    /* it is necesary to update computation of soundspeed if one uses
       alphaviscosity in FViscosity. It is not necesary in locally
       isothermal runs since cs is constant.  It is computed here for
       the needs of ConditionCFL. */
  }
  if (IsDisk == YES) {
    CommunicateBoundaries (Rho, Vrad, Vtheta, Energy, Label);
    if (SloppyCFL == YES)
      gastimestepcfl = ConditionCFL (Vrad, Vtheta,Rho,Energy, DT-dtemp);
  }
  MPI_Allreduce (&gastimestepcfl, &GasTimeStepsCFL, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  dt = DT / (real)GasTimeStepsCFL;
  while (dtemp < 0.999999999*DT) {
    MassTaper = PhysicalTime/(MASSTAPER*2.0*M_PI);
    MassTaper = (MassTaper > 1.0 ? 1.0 : pow(sin(MassTaper*M_PI/2.0),2.0));
    if (IsDisk == YES) {
      CommunicateBoundaries (Rho, Vrad, Vtheta, Energy, Label);
      if (SloppyCFL == NO) {
	gastimestepcfl = 1;
	gastimestepcfl = ConditionCFL (Vrad, Vtheta, Rho, Energy, DT-dtemp);
	MPI_Allreduce (&gastimestepcfl, &GasTimeStepsCFL, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
	dt = (DT-dtemp)/(real)GasTimeStepsCFL;
      }
      AccreteOntoPlanets (Rho, Vrad, Vtheta, dt, sys);
    }
    dtemp += dt;
    DiskOnPrimaryAcceleration.x = 0.0;
    DiskOnPrimaryAcceleration.y = 0.0;
    if (Corotating == YES) GetPsysInfo (sys, MARK);
    if (IsDisk == YES) {
      /* Indirect term star's potential computed here */
      DiskOnPrimaryAcceleration   = ComputeAccel (force, Rho, 0.0, 0.0, 0.0, 0.0);
      /* Gravitational potential from star and planet(s) is computed
	 and stored here*/
      FillForcesArrays (sys, Rho, Energy);
      /* Planets' velocities are updated here from gravitationnal
	 interaction with disk */
      AdvanceSystemFromDisk (force, Rho, Energy, sys, dt);
    }
    /* Planets' positions and velocities are updated from
       gravitational interaction with star and other planets */
    AdvanceSystemRK5 (sys, dt);
    /* Below we correct vtheta, planet's position and velocities if we
       work in a frame non-centered on the star */
    if (Corotating == YES) {
      OmegaNew = GetPsysInfo(sys, GET) / dt;
      domega = OmegaNew-OmegaFrame;
      if (IsDisk == YES)
	CorrectVtheta (Vtheta, domega);
      OmegaFrame = OmegaNew;
    }
    RotatePsys (sys, OmegaFrame*dt);
    /* Now we update gas */
    if (IsDisk == YES) {
      ApplyBoundaryCondition (Vrad, Vtheta, Rho, Energy, dt);
      Crashed = DetectCrash (Rho);    /* test for negative density values */
      Crashed = DetectCrash (Energy);  /* test for negative energy values */
      if (Crashed == YES && 0) {
	if (AlreadyCrashed == 0) {
          WriteDiskPolar (Rho, 9999);    /* We write the HD arrays */
          WriteDiskPolar (Vrad, 9999);   /* in order to keep a track */
          WriteDiskPolar (Vtheta, 9999); /* of what happened */
          WriteDiskPolar (Temperature, 9999);
	  timeCRASH=PhysicalTime;   /* if it appears to be the first crash */
	  fprintf (stdout,"\nCrash! at time %.12g\n", timeCRASH);
	  fflush (stdout);
	}
	AlreadyCrashed++;
	masterprint ("c");
      } else {
	masterprint (".");
      }
      fflush (stdout);
      if (ZMPlus)
	compute_anisotropic_pressurecoeff (sys);
      ComputePressureField (Rho, Energy);
      SubStep1 (Vrad, Vtheta, Rho, dt);
      SubStep2 (Rho, Energy, dt);
      ActualiseGas (Vrad, VradNew);
      ActualiseGas (Vtheta, VthetaNew);
      ApplyBoundaryCondition (Vrad, Vtheta, Rho, Energy, dt);//1st
      if (Adiabatic) {
	ComputeViscousTerms (Vrad, Vtheta, Rho);
	SubStep3 (Rho, dt);
	ActualiseGas (Energy, EnergyNew);
      }
      Transport (Rho, Vrad, Vtheta, Energy, Label, dt);
      ApplyBoundaryCondition (Vrad, Vtheta, Rho, Energy, dt);//2nd
      ComputeTemperatureField (Rho, Energy);
      ApplyBoundaryCondition (Vrad, Vtheta, Rho, Energy, dt);//3rd

      mdcp = CircumPlanetaryMass (Rho, sys);
      exces_mdcp = mdcp - mdcp0;
    }
    PhysicalTime += dt;
  }
  masterprint ("\n");
}

void Coolingf(Energy, Rho, dt)
     PolarGrid *Energy, *Rho;
     real dt;
{
  int i, j, l, nr, ns;
  real *energ, *qcoola, *diag3, *rho, newdt, subdt, lsun;
  double newdt2;
  nr = Energy->Nrad;
  ns = Energy->Nsec;
  energ = Energy->Field;
  rho = Rho->Field;
  qcoola = Qcoola->Field;
  //diag3 = Diag3->Field;
  subdt = 0.;
  lsun = 3.839e33/EUNIT*TUNIT;

  if(SUBCOOL){
   while (subdt < dt){
    //newdt=SubTimeStep(Energy, Rho);
    newdt=1e30;
    if(newdt>(dt-subdt)) newdt=dt-subdt;
    subdt += newdt;
    qcoola = Qcoola->Field;
    for (i = 0; i < nr; i++) {
     for (j = 0; j < ns; j++) {
      l = j+i*ns;
      energ[l] -= qcoola[l]*newdt;
     }
    }
   }
  }else{
   for (i = 0; i < nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      energ[l] -= qcoola[l]*dt;
    }
   }
  }
}


real SubTimeStep(Energy, Rho)
     PolarGrid *Energy,*Rho;
{
    int nr, ns, i, j, l;
    real *energ, *qcoola, *rho, dt, newdt, tem;
    nr = Energy->Nrad;
    ns = Energy->Nsec;
    energ = Energy->Field;
    qcoola = Qcoola->Field;
    rho = Rho->Field;
    newdt = 1.e30;
    CoolInit(Energy, Rho);
    for (i = 0; i < nr; i++) {
     for (j = 0; j < ns; j++) {
      l = j+i*ns;
      dt = fabs(0.3*energ[l]/qcoola[l]);
      if(dt < newdt) newdt = dt;
     }
    }
    return newdt;
}



void SubStep1 (Vrad, Vtheta, Rho, dt)
     PolarGrid *Vrad, *Vtheta, *Rho;
     real dt;
{
  int i, j, l, lim, ljm, ljp, nr, ns;
  boolean selfgravityupdate;
  extern boolean Evanescent;
  real *vrad, *vtheta, *rho;
  real *Pot, *Press;
  real *vradint, *vthetaint;
  real gradp, gradphi, vt2, dxtheta;
  real invdxtheta;
  real supp_torque=0.0;		/* for imposed disk drift */
  nr = Vrad->Nrad;
  ns = Vrad->Nsec;
  rho = Rho->Field;
  vrad = Vrad->Field;
  vtheta = Vtheta->Field;
  vradint = VradInt->Field;
  vthetaint = VthetaInt->Field;
  Pot = Potential->Field;
  Press = Pressure->Field;
  /* In this substep we take into account  */
  /* the source part of Euler equations  */
  /* We evolve velocities with pressure gradients */
  /* gravitational forces and curvature terms */
#pragma omp parallel private(j,l,lim,ljm,ljp,dxtheta,vt2,gradp,gradphi,invdxtheta,supp_torque)
  {
#pragma omp for nowait
    for (i = 1; i < nr; i++) {
      for (j = 0; j < ns; j++) {
	l = j+i*ns;
	lim = l-ns;
	ljp = l+1;
	if (j == ns-1) ljp = i*ns;
	gradp = (Press[l]-Press[lim])*2.0/(rho[l]+rho[lim])*InvDiffRmed[i];
	gradphi = (Pot[l]-Pot[lim])*InvDiffRmed[i];
	vt2 = vtheta[l]+vtheta[ljp]+vtheta[lim]+vtheta[ljp-ns];
	vt2 = vt2/4.0+Rinf[i]*OmegaFrame;
	vt2 = vt2*vt2;
	vradint[l] = vrad[l]+dt*(-gradp-gradphi+vt2*InvRinf[i]);
      }
    }
#pragma omp for
    for (i = 0; i < nr; i++) {
      supp_torque = IMPOSEDDISKDRIFT*0.5*pow(Rmed[i],-2.5+SIGMASLOPE);
      dxtheta = 2.0*PI/(real)ns*Rmed[i];
      invdxtheta = 1.0/dxtheta;
      for (j = 0; j < ns; j++) {
	l = j+i*ns;
	ljm = l-1;
	if (j == 0) ljm = i*ns+ns-1;
	gradp = (Press[l]-Press[ljm])*2.0/(rho[l]+rho[ljm])*invdxtheta;
	if ( ZMPlus )
	  gradp *= SG_aniso_coeff;
	gradphi = (Pot[l]-Pot[ljm])*invdxtheta;
	vthetaint[l] = vtheta[l]-dt*(gradp+gradphi);
	vthetaint[l] += dt*supp_torque;
      }
    }
  }
#ifdef SelfGravity
//  if ( SelfGravity ) {
    selfgravityupdate = YES;
    compute_selfgravity(Rho, VradInt, VthetaInt, dt, selfgravityupdate);
  //}
#endif
  ComputeViscousTerms (VradInt, VthetaInt, Rho);
  UpdateVelocitiesWithViscosity (VradInt, VthetaInt, Rho, dt);
  if ( !Evanescent )
    ApplySubKeplerianBoundary (VthetaInt);
    //Mdisk = GasTotalMass(Rho);

}

void SubStep2 (Rho, Energy, dt)
     PolarGrid *Rho, *Energy;
     real dt;
{
  int i, j, l, lim, lip, ljm, ljp, nr, ns;
  real *vrad, *vtheta, *rho, *energy;
  real *vradnew, *vthetanew, *qt, *qr, *energyint;
  real dxtheta, invdxtheta;
  real dv;
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  rho = Rho->Field;
  vrad = VradInt->Field;
  vtheta = VthetaInt->Field;
  qr = RhoInt->Field;
  vradnew = VradNew->Field;
  vthetanew = VthetaNew->Field;
  qt = TemperInt->Field;
  energy = Energy->Field;
  energyint = EnergyInt->Field;
#pragma omp parallel for private(j,dxtheta,l,lim,lip,ljm,ljp,dv)
  for (i = 0; i < nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      lip = l+ns;
      ljp = l+1;
      if (j == ns-1) ljp = i*ns;
      dv = vrad[lip]-vrad[l];
      if (dv < 0.0)
        qr[l] = CVNR*CVNR*rho[l]*dv*dv;
      else 
        qr[l] = 0.0; 
      dv = vtheta[ljp]-vtheta[l];
      if (dv < 0.0)
        qt[l] = CVNR*CVNR*rho[l]*dv*dv;
      else
	qt[l] = 0.0;
	//if (i==0 && j==0 &&) printf("vt=%e\t vr=%e\n",vtheta[l],vrad[l]);
    }
  }
#pragma omp parallel private(l,lim,lip,ljm,ljp,j,dxtheta,invdxtheta)
  {
#pragma omp for nowait
    for (i = 1; i < nr; i++) {
      for (j = 0; j < ns; j++) {
	l = j+i*ns;
	lim = l-ns;
	vradnew[l] = vrad[l]-dt*2.0/(rho[l]+rho[lim])*(qr[l]-qr[lim])*InvDiffRmed[i];
      }
    }
#pragma omp for
    for (i = 0; i < nr; i++) {
      dxtheta = 2.0*PI/(real)ns*Rmed[i];
      invdxtheta = 1.0/dxtheta;
      for (j = 0; j < ns; j++) {
	l = j+i*ns;
	ljm = l-1;
	if (j == 0) ljm = i*ns+ns-1;
	vthetanew[l] = vtheta[l]-dt*2.0/(rho[l]+rho[ljm])*(qt[l]-qt[ljm])*invdxtheta;
      }
    }
    /* If gas disk is adiabatic, we add artificial viscosity as a source */
    /* term for advection of thermal energy polargrid */
    if (Adiabatic) {
#pragma omp for nowait
      for (i = 0; i < nr; i++) {
	dxtheta = 2.0*PI/(real)ns*Rmed[i];
	invdxtheta = 1.0/dxtheta;
	for (j = 0; j < ns; j++) {
	  l = j+i*ns;
	  lip = l+ns;
	  ljp = l+1;
	  if (j == ns-1) ljp = i*ns;
	  energyint[l] = energy[l] -				\
	    dt*qr[l]*(vrad[lip]-vrad[l])*InvDiffRsup[i] -	\
	    dt*qt[l]*(vtheta[ljp]-vtheta[l])*invdxtheta;
	}
      }
    }
  }
}
	       
void SubStep3 (Rho, dt)
     PolarGrid *Rho;
     real dt;
{
  extern boolean Cooling;
  int i, j, l, nr, ns, k;
  int lip, li2p;
  real r, rip, ri2p, qpip, qpi2p, ff, tau, taua,qp1,qp2;
  real *dens, *pres, *energy, *energynew, rand;
  real *divergence, *Trr, *Trp, *Tpp, *qplus, *opa, *qcoola;//, *diag0, *diag1, *diag2, *diag3, *diag4;
  real viscosity, energypred, num, den, frac, CoolingTime;
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  dens = Rho->Field;
  pres = Pressure->Field;
  energy = EnergyInt->Field;
  energynew = EnergyNew->Field;
  divergence = DivergenceVelocity->Field;
  qplus = Qplus->Field;
  qcoola = Qcoola->Field;

  Trr = TAURR->Field;
  Trp = TAURP->Field;
  Tpp = TAUPP->Field;

//  Trr2 = TAURR2->Field;
 // Trp2 = TAURP2->Field;
 // Tpp2 = TAUPP2->Field;
  opa = Opacity->Field;

  /* In this substep we take into account  */
  /* the source part of energy equation  */
  /* We evolve internal energy with  */
  /* compression/dilatation and heating terms */
#pragma omp parallel private(j,l,viscosity,energypred)
  {
#pragma omp for
    /* We calculate the heating source term Qplus from i=1 */
    for (i = 1; i < nr; i++) { /* Trp defined from i=1 */
      viscosity = FViscosity (Rmed[i]); /* Global_Bufarray holds cs */

      k=0;
      while (GlobalRmed[k] < Rmed[i]) k++;

      for (j = 0; j < ns; j++) {
	l = j+i*ns;
         //frac = min2(SIGMAA,dens[l])/dens[l];
	frac=1.0;
	if (viscosity != 0.0) {
          tau = 0.5*opa[l]*dens[l];
          //taua = 0.5*opa[l]*min2(dens[l],SIGMAA);
	  taua= 0.5*opa[l]*dens[l];
	  qplus[l] = 0.5*frac*(taua/tau)/viscosity/dens[l]*(1.0*Trr[l]*Trr[l] +	\
					   2.0* Trp[l]*Trp[l] +	\
					    1.0*Tpp[l]*Tpp[l] );
	qp1=qplus[l];
	  //qplus[l] += 0.5*(1-frac)/viscosity/dens[l]*(Trr[l]*Trr[l] +  \
                                            Trp[l]*Trp[l] +     \
                                            Tpp[l]*Tpp[l] );
           //diag1[l] = qplus[l];
	qp2=(2.0/9.0)*frac*(taua/tau)*viscosity*dens[l]*divergence[l]*divergence[l];

           qplus[l] += (2.0/9.0)*frac*(taua/tau)*viscosity*dens[l]*divergence[l]*divergence[l];
           //qplus[l] += (2.0/9.0)*(1.-frac)*viscosity*dens[l]*divergence[l]*divergence[l];
          // diag2[l] = qplus[l] - diag1[l];
	//if (i==10 && j==10  && 1) printf("Trr=%e \t Trp= %e\t Tpp= %e\t div= %e\t r=%f\t qp1=%e\t qp2=%e\n ",Trr[l],Trp[l],Tpp[l],divergence[l],Rmed[i],qp1,qp2);
	qplus[l]=qplus[l];
	}
	else
	  qplus[l] = 0.0;
      }
//printf("%f\t%1.20f\n",Rmed[i],qplus[l]);
    }

    /* We calculate the heating source term Qplus for i=0 */
    i = 0;
    r    = Rmed[i];
    rip  = Rmed[i+1];
    ri2p = Rmed[i+2];
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      lip = l+ns;
      li2p = lip+ns;
      qpip = qplus[lip];   /* qplus(i=1,j) */
      qpi2p = qplus[li2p]; /* qplus(i=2,j) */
      if (viscosity != 0.0) {
	/* power-law extrapolation */
	qplus[l] = qpip*exp( log(qpip/qpi2p) * log(r/rip) / log(rip/ri2p) );
       // diag1[l] = qplus[l];
      }
      else
	qplus[l] = 0.0;
    }
    /* Now we can update energy with source terms from i=0 */
    for (i = 0; i < nr; i++) {
      for (j = 0; j < ns; j++) {
	l = j+i*ns;
	if (!Cooling) {
	  num = dt*qplus[l] + energy[l];
	  den = 1.0+(ADIABATICINDEX-1.0)*dt*divergence[l];
	  energynew[l] = num/den;
	}
	if (Cooling) {
          if (COOLMETHOD==1){
          if(PhysicalTime<30.) {
              //  ff=29*(30.-PhysicalTime)/30.+1. ;
		ff=1.;
          }else{
                ff=1.;
          }
          num = ff*EnergyMed[i]*dt*dens[l]/SigmaMed[i] + CoolingTimeMed[i]*energy[l] + dt*CoolingTimeMed[i]*(qplus[l]);
          den = dt + CoolingTimeMed[i] + (ADIABATICINDEX-1.0)*dt*CoolingTimeMed[i]*divergence[l];
          //printf("%1.12f\t%1.12f\n",energy[l],qcoola[l]);
	  energynew[l] = num/den;
         }
         if (COOLMETHOD==2){
	 //CoolInit(Energy, Rho);
          //num = dt*qplus[l] + energy[l];
	  //if ((energy[l] < 0.0 || 1) && i==0 && j==0) printf("%g\t%g\t%g\t qplus=%g\t qcoola=%g\t %g\t dt=%g\t %d\t%d\n",dt*qplus[l],energy[l],dt*qcoola[l],qplus[l],qcoola[l],Trp[l],dt,i,j);
          num = dt*qplus[l] - dt*qcoola[l] + energy[l];
          den = 1.0+(ADIABATICINDEX-1.0)*dt*divergence[l];
          energynew[l] = num/den;
     
   

    // diag3[l] = (ADIABATICINDEX-1.0)*energy[l]*divergence[l];

          /*if(ORBITCOOLINNER){
            if(Rmed[i]<2.5*RMIN){
              CoolingTime=CoolingTimeMed[i]*((Rmed[i]/RMIN-1.)*10.+1.);
              num = EnergyMed[i]*dt*dens[l]/SigmaMed[i] + CoolingTime*energy[l] + dt*CoolingTime*(qplus[l]);
              den = dt + CoolingTime + (ADIABATICINDEX-1.0)*dt*CoolingTime*divergence[l];
              energynew[l] = num/den;
            }
          }*/
         }
	}
      }
    }
  }
}
 

void CoolInit(Energy, Rho)
     PolarGrid *Energy,*Rho;
{
    int nr, ns, i, j, l, lim, ljm, lip, ljp, istart, iend;
    real *energ, *opa, *rho,  *qcoola, *temp,*opar; 
    real hei, heia, rhoc, rhoccgs, tccgs, kappacgs, tau, taua, Text4, Lstarcgs, Lacccgs, alpha, dens_temp, opac0;
    real lsun, dxtheta, invdxtheta, dphi, invdphi;
    nr = Energy->Nrad;
    ns = Energy->Nsec;
    energ = Energy->Field;
    opa = Opacity->Field;
   // opar = Opacityr->Field;
    qcoola = Qcoola->Field;
    rho = Rho->Field;
    //tdr = Tdr->Field;
    //fluxr = Fluxr->Field;
    //fluxt = Fluxt->Field;
    temp = Temperature->Field;
    //diag1 = Diag1->Field;
    //diag2 = Diag2->Field;
    //diag3 = Diag3->Field;
    //diag4 = Diag4->Field;
    lsun = 3.839e33/EUNIT*TUNIT;
    dphi = 2.0*M_PI/(real)ns;
    invdphi = 1.0/dphi;
    Lstarcgs = LSUN*pow(10.,0.2)*pow(MC*(MUNIT/MSUN),1.74);
    Lstar = Lstarcgs*TUNIT/EUNIT;
   // Lacccgs = 0.5*xgconstcgs*(MC*MUNIT)*(fabs(AccRate)*MUNIT/TUNIT)/RSUN;
   // Lacc = Lacccgs*TUNIT/EUNIT;

   /* if (RADIALCOOL){
      istart=1;
      iend=nr-1;
    }else{
      istart=0;
      iend=nr;
    }*/
      istart=0;
      iend=nr;
 
    for (i = 0; i < nr; i++) {
	opac0=OpacMed[i];
      for (j = 0; j < ns; j++) {
	  l = j+i*ns;
if (Restart && 0) {
        energ[l]=EnergyMed[i];
        }
        hei=sqrt(ADIABATICINDEX*(ADIABATICINDEX-1.0)*(energ[l]/rho[l]))*pow(Rmed[i],1.5)/sqrt(G*MC);
        rhoc=rho[l]/sqrt(2.*M_PI)/hei;
        rhoccgs = rhoc*MUNIT/LUNIT/LUNIT/LUNIT;
        temp[l] = MU/R*(ADIABATICINDEX-1.0)*energ[l]/rho[l];
        //tccgs = temp[l]*TEUNIT;
        tccgs=temp[l]*LUNIT*LUNIT/TUNIT/TUNIT/8.31446e7; 
	//opa[l]=opac(max2(1.0e-4/2./hei,rhoccgs), (tccgs));
 	
	if (TimeStep>1) 
	{
	opa[l]=opac((rhoccgs), (tccgs));
}
	else{
	opa[l]=6.0;
	}
       //opa[l]=10.0;
	opa[l]=opa[l]/LUNIT/LUNIT*MUNIT;
   //     opar[l]=opa[l]*rhoc;
        //tdr[l]=-16./3.*sigma*temp[l]*temp[l]*temp[l]*temp[l]*2.*hei;
      }
    }
  /*  if(RADIALCOOL){
     for (i = 1; i < nr; i++) {
      for (j = 0; j < ns; j++) {
        l = j+i*ns;
        lim = l-ns;
        fluxr[l] = 0.5*(tdr[l]+tdr[lim])*(temp[l]-temp[lim])/(1./InvDiffRmed[i]*0.5*(opar[l]+opar[lim])+InvDiffRmed[i]/0.5/(opar[l]+opar[lim]));
      }
     }
     for (i = 0; i < nr; i++) {
      for (j = 0; j < ns; j++) {
        l = j+i*ns;
        ljm = l-1;
        if (j == 0) ljm = i*ns+ns-1;
        dxtheta = 2.0*PI/(real)ns*Rmed[i];
        invdxtheta = 1.0/dxtheta;
        fluxt[l] = 0.5*(tdr[l]+tdr[ljm])*(temp[l]-temp[ljm])/(1./invdxtheta*0.5*(opar[l]+opar[ljm])+invdxtheta/0.5/(opar[l]+opar[ljm]));
      }
     }
    }
*/
    for (i = istart; i < iend; i++) {
      for (j = 0; j < ns; j++) {
        l = j+i*ns;
        lip=l+ns;
        ljp=l+1;
        if (j == ns-1) ljp = i*ns;

        tau = 0.5*rho[l]*opa[l];
  //      diag4[l] = tau;

        //if (rho[l] <= SIGMAA) alpha = ALPHAMRI;
        //else alpha = (ALPHAMRI*SIGMAA+ALPHARD*(rho[l]-SIGMAA))/rho[l];

	//alpha=ALPHAVISCOSITY;
        //taueff = 3.*tau/8.*max2(1.0e-10,(1.+min2(SIGMAA,rho[l])/rho[l]-ALPHAMRI*min2(SIGMAA,rho[l])/alpha/rho[l]))+sqrt(3.)/4.+1./3./tau;

        //taua = 0.5*(rho[l],SIGMAA)*opa[l];
        taua=0.5*rho[l]*opa[l];
        //Text4 = 0.1*(Lstar+Lacc)/4./PI/Rmed[i]/Rmed[i]/sigma+pow(TENV,4.);
        Text4 = 0.1*Lstar/4./PI/Rmed[i]/Rmed[i]/sigma;

	qcoola[l] = 16./3.*sigma*(temp[l]*temp[l]*temp[l]*temp[l]-Text4)*tau/(1.+tau*tau);
        //printf("%1.20f\n",qcoola[l]);
	/* if(RADIALCOOL){
        qcoola[l] += (fluxr[lip]*Rsup[i]-fluxr[l]*Rinf[i])*InvDiffRsup[i]*InvRmed[i];
        qcoola[l] += (fluxt[ljp]-fluxt[l])*invdphi*InvRmed[i];
        }*/
      }
    }

}


	   
int ConditionCFL (Vrad, Vtheta, Rho, Energy, deltaT)
     PolarGrid *Vrad, *Vtheta, *Rho, *Energy;
     real deltaT;
{
  static real Vresidual[MAX1D], Vmoy[MAX1D];
  int i, j, l, ns, nr, lip, ljp;
  real invdt1, invdt2, invdt3, invdt4, cs, newdt, dt;
  int ideb, jdeb;
  real itdbg1, itdbg2, itdbg3, itdbg4,itdbg5, mdtdbg; /* debugging variables */
  real *vt, *vr, *rho, *energ, *qcoola, dxrad, dxtheta, dvr, dvt, viscr, visct, tem, tema;
  real *soundspeed;
  soundspeed = SoundSpeed->Field;
  ns = Vtheta->Nsec;
  nr = Vtheta->Nrad;
  vt = Vtheta->Field;
  vr = Vrad->Field;
  rho = Rho->Field;
  qcoola = Qcoola->Field;
  energ = Energy->Field;
  newdt = 1e30;
  for (i = 0; i < nr; i++) {
    dxrad = Rsup[i]-Rinf[i];
    dxtheta=Rmed[i]*2.0*PI/(real)ns;
    Vmoy[i] = 0.0;
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      Vmoy[i] += vt[l];
    }
    Vmoy[i] /= (real)ns;
  }
  for (i = One_or_active; i < Max_or_active; i++) {
    dxrad = Rsup[i]-Rinf[i];
    dxtheta=Rmed[i]*2.0*PI/(real)ns;
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      if (FastTransport == YES)
	Vresidual[j] = vt[l]-Vmoy[i];  /* FARGO algorithm */
      else
	Vresidual[j] = vt[l];	       /* Standard algorithm */
    }
    Vresidual[ns]=Vresidual[0];
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      lip = l+ns;
      ljp = l+1;
      if (j == ns-1) ljp = i*ns;
      cs = soundspeed[l];
      invdt1 = cs/(min2(dxrad,dxtheta));
      invdt2 = fabs(vr[l])/dxrad;
      invdt3 = fabs(Vresidual[j])/dxtheta;
      dvr = vr[lip]-vr[l];
      dvt = vt[ljp]-vt[l];
      if (dvr >= 0.0) dvr = 1e-10;
      else dvr = -dvr;
      if (dvt >= 0.0) dvt = 1e-10;
      else dvt = -dvt;
      invdt4 = max2(dvr/dxrad,dvt/dxtheta);
      invdt4*= 4.0*CVNR*CVNR;
      dt = CFLSECURITY/sqrt(invdt1*invdt1+invdt2*invdt2+invdt3*invdt3+invdt4*invdt4);
      if (dt < newdt) {
	newdt = dt;
	if (debug == YES) {
	  ideb = i;
	  jdeb = j;
	  itdbg1 = 1.0/invdt1; itdbg2=1.0/invdt2; itdbg3=1.0/invdt3; itdbg4=1.0/invdt4;
	  mdtdbg = newdt;
	  viscr = dxrad/dvr/4.0/CVNR/CVNR;     
	  visct = dxtheta/dvt/4.0/CVNR/CVNR;
	}
      }  
   }
  }
  for (i = Zero_or_active; i < MaxMO_or_active; i++) {
    dt = 2.0*PI*CFLSECURITY/(real)NSEC/fabs(Vmoy[i]*InvRmed[i]-Vmoy[i+1]*InvRmed[i+1]);
    if (dt < newdt) newdt = dt;
  }

  if (COOLMETHOD==2){
    CoolInit(Energy, Rho);
    for (i = Zero_or_active; i < Max_or_active; i++) {
      for (j = 0; j < ns; j++) {
        l = j+i*ns;
        if(SUBCOOL){
          dt= pow((newdt/fabs(0.3*energ[l]/qcoola[l])),1.)*fabs(0.3*energ[l]/qcoola[l]);
        }else{
          dt= fabs(0.3*energ[l]/qcoola[l]);
	   // dt = 1e20;
        }
        if (dt < newdt) {
          newdt = dt;
//        if (debug == YES) {
            ideb = i;
            jdeb = j;
            itdbg5 = dt;
//        }
        }
      }
    }
  }

  if (deltaT < newdt) newdt = deltaT;
  if (debug == YES) {
    printf ("Timestep control information for CPU %d: \n", CPU_Rank);
    printf ("Most restrictive cell at i=%d and j=%d\n", ideb, jdeb);
    printf ("located at radius Rmed         : %g\n", Rmed[ideb]);
    printf ("Sound speed limit              : %g\n", itdbg1);
    printf ("Radial motion limit            : %g\n", itdbg2);
    printf ("Residual circular motion limit : %g\n", itdbg3);
    printf ("Viscosity limit                : %g\n", itdbg4);
    printf ("radiative cooling              : %g\n", itdbg5);
    printf ("   Arise from r with limit     : %g\n", viscr);
    printf ("   and from theta with limit   : %g\n", visct);
    printf ("Limit time step for this cell  : %g\n", mdtdbg);
    printf ("Limit time step adopted        : %g\n", newdt);
    if (newdt < mdtdbg) {
      printf ("Discrepancy arise either from shear.\n");
      printf ("or from the imposed DT interval.\n");
    }
  fflush (stdout);

  }
  return (int)(ceil(deltaT/newdt));
}

void ComputeSoundSpeed (Rho, Energy)
     PolarGrid *Rho;
     PolarGrid *Energy;
{
  int i, j, l, nr, ns, ii;
  real *dens, *energ, *cs;
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  dens = Rho->Field;
  energ = Energy->Field;
  cs = SoundSpeed->Field;
  for ( i = 0; i < nr; i++ ) {
    for ( j = 0; j < ns; j++ ) {
      l = i*ns + j;
      if (!Adiabatic)
	cs[l] = AspectRatio(Rmed[i])*sqrt(G*MC/Rmed[i])*pow(Rmed[i], FLARINGINDEX);
      else
	cs[l] = sqrt( ADIABATICINDEX*(ADIABATICINDEX-1.0)*energ[l]/dens[l] );
    }
  }
}

void ComputePressureField (Rho, Energy)
     PolarGrid *Rho;
     PolarGrid *Energy;
{
  int i, j, l, nr, ns, ii;
  real *dens, *pres, *energ, *cs;
  real peq, dp;
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  dens = Rho->Field;
  pres = Pressure->Field;
  energ = Energy->Field;
  cs = SoundSpeed->Field;
  for ( i = 0; i < nr; i++ ) {
    for ( j = 0; j < ns; j++ ) {
      l = i*ns + j;
      if (!Adiabatic) {
	pres[l] = dens[l]*cs[l]*cs[l]; /* since SoundSpeed is not updated */
                                       /* from initialization, cs remains */ 
                                       /* axisymmetric */
      }
      else
	pres[l] = (ADIABATICINDEX-1.0)*energ[l];
    }
  }
}

void ComputeTemperatureField (Rho, Energy)
     PolarGrid *Rho, *Energy;
{
  int i, j, l, nr, ns;
  real *dens, *pres, *energ, *temp;
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  dens = Rho->Field;
  pres = Pressure->Field;
  energ = Energy->Field;
  temp = Temperature->Field;
  for ( i = 0; i < nr; i++ ) {
    for ( j = 0; j < ns; j++ ) {
      l = i*ns + j;
      if (!Adiabatic)
	temp[l] = MU/R* pres[l]/dens[l];
      else{
	temp[l] = MU/R*(ADIABATICINDEX-1.0)*energ[l]/dens[l];
	//printf("temp=%f\n",temp[l]);
    }
	}
  }
}

real CircumPlanetaryMass (Rho, sys)
     PolarGrid *Rho;
     PlanetarySystem *sys;
{
  int i, j, l, ns;
  real xpl, ypl;
  real dist, mdcplocal, mdcptotal;
  real *dens, *abs, *ord;
  ns = Rho->Nsec;
  dens = Rho->Field;
  abs = CellAbscissa->Field;
  ord = CellOrdinate->Field;
  xpl = sys->x[0];
  ypl = sys->y[0];
  mdcplocal = 0.0;
  mdcptotal = 0.0;
  if (FakeSequential && (CPU_Rank > 0)) 
    MPI_Recv (&mdcplocal, 1, MPI_DOUBLE, CPU_Rank-1, 0, MPI_COMM_WORLD, &stat);
  for ( i = Zero_or_active; i < Max_or_active; i++ ) {
    for ( j = 0; j < ns; j++ ) {
      l = i*ns + j;
      dist = sqrt ( (abs[l]-xpl)*(abs[l]-xpl) +		\
		    (ord[l]-ypl)*(ord[l]-ypl) );
      if ( dist < HillRadius ) {
	mdcplocal += Surf[i] * dens[l];
      }
    }
  }
  if (FakeSequential) {
     if (CPU_Rank < CPU_Number-1)
       MPI_Send (&mdcplocal, 1, MPI_DOUBLE, CPU_Rank+1, 0, MPI_COMM_WORLD);
  }
  else
    MPI_Allreduce (&mdcplocal, &mdcptotal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  if (FakeSequential) {
    MPI_Bcast (&mdcplocal, 1, MPI_DOUBLE, CPU_Number-1, MPI_COMM_WORLD);
    mdcptotal = mdcplocal;
  }
  return mdcptotal;
}
