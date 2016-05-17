
/* generate opacity table from power-law fits */

/*

$Id: opac.c,v 3.2 2007/02/23 15:35:01 gammie Exp $

*/
#include "types.h"
#include <stdio.h>
#include <math.h>
real opac(real rho, real T)
{
	real xlop, xlp, xlt, rosstabd,pre;
	pre=rho*T*8.314472*1.e7/2.4;
	if (pre < 0. || T < 0.) {
		fprintf(stderr, "error: pre or T negative\n");
		fprintf(stderr, "rho: %g pre: %g T: %g\n", rho, pre, T);
		prs_exit (1);
	}
	if (pre == 0 || T == 0)
		return (1.);
	xlp = log10(pre);
	xlt = log10(T);


        if (xlt<3.+0.03*(xlp+4.)){
          xlop=-1.27692+0.73846*xlt;
	}
        else if (xlt<3.08+0.028084*(xlp+4)){
          xlop=129.88071-42.98075*xlt+(142.996475-129.88071)*0.1*(xlp+4);
	}
        else if (xlt<3.28+xlp/4.*0.12){
          xlop=-15.0125+4.0625*xlt;
	}
        else if (xlt<3.41+0.03328*xlp/4.){
          xlop=58.9294-18.4808*xlt+(61.6346-58.9294)*xlp/4.;
	}
        else if (xlt<3.76+(xlp-4)/2.*0.03){
          xlop=-12.002+2.90477*xlt+(xlp-4)/4.*(13.9953-12.002);
	}
        else if (xlt<4.07+(xlp-4)/2.*0.08){
          xlop=-39.4077+10.1935*xlt+(xlp-4)/2.*(40.1719-39.4077);
	}
        else if (xlt<5.3715+(xlp-6)/2.*0.5594){
          xlop=17.5935-3.3647*xlt+(xlp-6)/2.*(17.5935-15.7376);
	}
        else{
        xlop=-0.48;
        }
        if (xlop<3.586*xlt-16.85&&xlt<4.){xlop=3.586*xlt-16.85;}
        if (xlt<2.9){xlop=-1.27692+0.73846*xlt;}
        rosstabd=pow(10.,xlop);
        return(rosstabd);

}

