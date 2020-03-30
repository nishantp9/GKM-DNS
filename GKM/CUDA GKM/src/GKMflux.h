/****************************************************************************************
 * Header for Flux calculation using the Gas-Kinetic scheme										 *
 * -------------------------------------------------------------------------------------*
 ****************************************************************************************/

#ifndef GKMFLUXHEADERDEF
#define GKMFLUXHEADERDEF

#include "param.h"
void flux(ptype W0[5], ptype DWx[5], ptype DWy[5],
			 ptype DWz[5],ptype mu, ptype dt, ptype dx, ptype F[5]);

void slopesolver(ptype b[5], ptype U[3], ptype lam, ptype a[5]);

void mult(ptype M[5][2], int k, ptype C, ptype B[5]);

void Mult(ptype M[5][5], ptype C, ptype A[5], ptype B[5]);

void c2p(ptype W[5], ptype &den, ptype &Ux, ptype &Uy, ptype &Uz, ptype &P);

void add(ptype A[5], ptype B[5]);


#endif

/****************************************************************************************
 * -------------------------------------------------------------------------------------*/
