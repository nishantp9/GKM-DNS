#ifndef CPUUPDATEHEADERDEF
#define CPUUPDATEHEADERDEF

#include "param.h"
#include "FlowField.h"
#include <cmath>
#include <iostream>

class CPUevolve
{
private:
	ptype *Fx[5], *Fy[5], *Fz[5];

	ptype deni, ui, vi, wi, Ti, Pi, Vsqr;
	ptype WL[5] ,WR[5] ,W0[5];
	ptype WLN[5],WRN[5],W0N[5];
	
	ptype WLS[5],WRS[5],W0S[5];
	ptype WLF[5],WRF[5],W0F[5];
	ptype WLB[5],WRB[5],W0B[5];
		
	ptype DWx[5],DWy[5],DWz[5],F[5];	
	
	void derivsX(Macro *macro);
	void derivsY(Macro *macro);
	void derivsZ(Macro *macro);
	void update(Macro *macro);
	void c2p();
	void PrintF(ptype *F[], int l);

public:
	CPUevolve();
	~CPUevolve();		
	void evolve(Macro *macro);
	void PrintStats(Macro *macro);

};
#endif
