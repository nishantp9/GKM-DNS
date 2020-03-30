/***************************************************************
 * Creates all variables required for GKM DNS						*
 * Defines all necessary subroutines									*
 **************************************************************/
#ifndef FLOWFIELDHEADERDEF
#define FLOWFIELDHEADERDEF
 
#include <cmath>
#include "param.h"

class Macro
{
public:
	ptype t;
	ptype dx, dt, maxU, maxC, mu;
	ptype var, den0, urms0, mu0, c0, p0, Mt, Re, k0, T0, Vsqr;

	ptype *den, *p, *u, *v, *w, *c, *M;
	ptype *W[5];
		
	Macro();
	~Macro();

   void timestep();
   void Initialize_DHIT();
   void PeriodicBC();
   
   void SaveVelocityField();
   void SaveVelocityField2();
   void SaveVelocityField3();
   
   
   

private:
   void p2c3D();
   void c2p3D();
};

#endif
