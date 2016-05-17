#include "mp.h"

#ifdef SelfGravity
void compute_selfgravity (Rho, RadialVelocity, AzimutalVelocity, DeltaT, SGUpdate)
     PolarGrid *Rho;
     PolarGrid *RadialVelocity;
     PolarGrid *AzimutalVelocity;
     real DeltaT;
     boolean SGUpdate;
{
  if ( SG_initcounter == 0 ) {
    init_sg ();
#ifndef  SGZeroMode
    //if ( !SGZeroMode )
      compute_fftkernel ();
#endif
  }
  /* Only mode m=0 of disk self-gravity is taken into account */
#ifdef SGZeroMode
//  if ( SGZeroMode )
    compute_SGZeroMode (Rho);
#endif  
/* Complete disk self-gravity treatment in this case */
#ifndef  SGZeroMode
  //if ( !SGZeroMode ) {
    compute_fftdensity (Rho);
    /* Here we compute radial and azimutal components of sg acceleration
       as a convolution product of reduced density and kernel arrays */
    compute_sgacc (Rho);
//  }
#endif
  if ( SGUpdate ) {
    /* Computes polar components of acceleration and
       updates values of vrad, vtheta at each step */
    update_sgvelocity (RadialVelocity, AzimutalVelocity, DeltaT);
  }
  SG_initcounter++;
}
#endif
