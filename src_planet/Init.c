#include "mp.h"

extern boolean Restart;
extern int     NbRestart;

void readopac()
{
FILE *fp;

char opac_input[256];
double value1;
int i,j,foo;
if (CPU_Rank > 0) MPI_Recv (&foo, 1, MPI_INT, CPU_Prev, 10, MPI_COMM_WORLD, &stat);

sprintf(opac_input, "%sopac.dat", OUTPUTDIR);
fp=fopen(opac_input,"r");
if (fp==NULL) return;
for (i=0;i<MAX1D;i++){
        /*if (CPU_Rank==0) {
                fscanf(fp,"%lf",&value1);
                for(j=0;j<CPU_Highest;j++){
                        MPI_Send (&value1,1, MPI_DOUBLE, j+1, 13, MPI_COMM_WORLD);
                }
        }
        if (CPU_Rank>0) {
                MPI_Recv (&value1, 1, MPI_DOUBLE, 0, 13, MPI_COMM_WORLD, &stat);
        }
        MPI_Barrier (MPI_COMM_WORLD);

        OpacMed[i]=value1;
*/
fscanf(fp,"%lf",&value1);
        OpacMed[i]=value1;


}
fclose(fp);
if (CPU_Rank < CPU_Highest) MPI_Send (&foo, 1, MPI_INT, CPU_Next, 10, MPI_COMM_WORLD);
        MPI_Barrier (MPI_COMM_WORLD);

}





void RefillVtheta (vtheta)
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


void ReadfromFile (array, fileprefix, filenumber)
     PolarGrid *array;
     char *fileprefix;
     int filenumber;
{
  int nr, ns, c, foo=0;
  real *field;
  char name[256];
  FILE *input;
  /* Simultaneous read access to the same file have been observed to
     give wrong results. */
  /* A sequential reading is imposed below. */
  /* If current CPU has a predecessor, wait for a message from him */
  if (CPU_Rank > 0) MPI_Recv (&foo, 1, MPI_INT, CPU_Prev, 10, MPI_COMM_WORLD, &stat);
  sprintf (name, "%s%s%d.dat", OUTPUTDIR, fileprefix, filenumber);
  input = fopen (name, "r");
  if (input == NULL) {
    fprintf (stderr, "WARNING ! Can't read %s. Restarting with t=0 settings.\n", name); 
    if (CPU_Rank < CPU_Highest) MPI_Send (&foo, 1, MPI_INT, CPU_Next, 10, MPI_COMM_WORLD);
    return;
  }
  field = array->Field;
  nr = array->Nrad;
  ns = array->Nsec;
  for (c = 0; c < IMIN; c++) {
    fread (field, sizeof(real), ns, input); 
    /* Can't read at once in order not to overflow 'field' */
  }
  fread (field, sizeof(real), nr*ns, input);
  fclose (input);
  /* Next CPU is waiting. Tell it to start now by sending the message
     that it expects */
  if (CPU_Rank < CPU_Highest) MPI_Send (&foo, 1, MPI_INT, CPU_Next, 10, MPI_COMM_WORLD);
  MPI_Barrier (MPI_COMM_WORLD);	/* previous CPUs do not touch anything
				   meanwhile */
}

void InitLabel (array, sys)
     PolarGrid *array;
     PlanetarySystem *sys;
{
  int nr, ns, i, j, l;
  real xp, yp, rp;
  real x, y, angle, distance, rhill;
  real *field;
  field = array->Field;
  nr = array->Nrad;
  ns = array->Nsec;
  xp = sys->x[0];
  yp = sys->y[0];
  rp = sqrt ( xp*xp + yp*yp );
  rhill = rp * pow( sys->mass[0]/3., 1./3 );
   /* Initialize label as you wish. In this example, label only takes
      into account fluid elements inside the planet's Hill Sphere */
  for (i = 0; i < nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j + i*ns;
      angle = (real)j/(real)ns*2.0*PI;
      x = Rmed[i] * cos(angle);
      y = Rmed[i] * sin(angle);
      distance = sqrt( (x - xp)*(x - xp) + (y - yp)*(y - yp) );
      if ( distance < rhill )
	field[l] = 1.0;
      else
	field[l] = 0.0;
    }
  }
}

void InitUnit()
{
  TUNIT=1./sqrt(6.67e-8/LUNIT/LUNIT/LUNIT*MUNIT);
  EUNIT=MUNIT*LUNIT/TUNIT*LUNIT/TUNIT;
  MOUNIT=MUNIT;
  TEUNIT=1./8.31451e7*EUNIT/MOUNIT;
}

void InitParam()
{
  real tenvcgs, rgascgs;
  real cscgs, OMEGA_break, OMEGA_C, ricgs, tcgs, ginf, rccgs, dmfalcgs, ttotcgs;

  MC=1.0;
  SIGMA0 = SIGMA0/MUNIT*LUNIT*LUNIT;
  rgascgs=8.3143e7;
  xgconstcgs=6.67259e-8;
  cscgs=sqrt(rgascgs*tenvcgs/MU);
//  OMEGA_break = 2.*sqrt(2.)*cscgs*cscgs*cscgs/xgconstcgs/(M_CLOUD*MSUN);
  ricgs=2.85e16;
  tcgs = PhysicalTime*TUNIT;
//  ginf = OMEGA_C*(0.975/2.*cscgs*tcgs+ricgs*M_CLOUD*(MSUN/MUNIT))*(0.975/2.*cscgs*tcgs+ricgs*M_CLOUD*(MSUN/MUNIT));
  ginf = OMEGA_C*(0.975/2.*cscgs*tcgs+ricgs)*(0.975/2.*cscgs*tcgs+ricgs);
  rccgs = ginf*ginf/xgconstcgs/(MC*MUNIT);
  sigma = 5.67051e-5/EUNIT*LUNIT*LUNIT*TEUNIT*TEUNIT*TEUNIT*TEUNIT*TUNIT;
  printf("%1.20f\n",EUNIT);
}



void Initialization (gas_density, gas_v_rad, gas_v_theta, gas_energy, gas_label, pla_sys)
     PolarGrid *gas_density, *gas_v_rad, *gas_v_theta, *gas_energy, *gas_label;
     PlanetarySystem *pla_sys;
{
  extern boolean Adiabatic;
  real *energ, *dens;
  int i, j, l, nr, ns;
  energ = gas_energy->Field;
  nr = gas_energy->Nrad;
  ns = gas_energy->Nsec;
  ReadPrevDim ();
  InitEuler (gas_v_rad, gas_v_theta, gas_density, gas_energy);
  InitLabel (gas_label, pla_sys);
  if (Restart == YES) {
    CheckRebin (NbRestart);
    MPI_Barrier (MPI_COMM_WORLD);
    /* Don't start reading before master has finished rebining... */
    /* It shouldn't be a problem though since a sequential read is */
    /* imposed in the ReadfromFile function below */
    mastererr ("Reading restart files...");
    fflush (stderr);
    ReadfromFile (gas_density, "gasdens", NbRestart);
    ReadfromFile (gas_v_rad, "gasvrad", NbRestart);
    ReadfromFile (gas_v_theta, "gasvtheta", NbRestart);
    if (Adiabatic) {
      ReadfromFile (gas_energy, "gasTemperature", NbRestart);
      /* ! gas_energy accounts for the gas temperature... */
      dens = gas_density->Field;
      for (i=0; i<nr; i++) {
	for (j=0; j<ns; j++) {
	  l = i*ns + j;
	  energ[l] = dens[l]*energ[l]/(ADIABATICINDEX-1.0);
	  if (energ[l] <=0.0) printf("ENERGY_RESTART: %g\t %d\t %d\n",energ[l],i,j);
	  /* this is e = dens*temp / (gamma-1) */
	}
      }
    }
    ReadfromFile (gas_label, "gaslabel", NbRestart);
    if (StoreSigma) {
	RefillSigma (gas_density);
	RefillVtheta (gas_v_theta);
	RefillVrad (gas_v_rad);}
	readopac();
    /* StoreEnergy = NO if Adiabatic = NO */
    if (1) RefillEnergy (gas_energy);
    fprintf (stderr, "done\n");
    fflush (stderr);
  }
  WriteDim (); 
}
