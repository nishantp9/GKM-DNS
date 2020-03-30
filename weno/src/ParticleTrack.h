/***************************************************************
 * Tracks the Particles
 **************************************************************/
#ifndef PARTICLETRACKHEADERDEF
#define PARTICLETRACKHEADERDEF

#include <fstream>
#include <iostream>
#include "param.h"
#include <bspline.h>

using namespace std;
class Particle
{
public:
	Particle();
	int infoSz, tagId; ptype info[8];
	bool blockChanged;
	ptype pos[3], DEN, U[3], P, gradD[3], gradU[3], gradV[3], gradW[3], gradP[3];
	
	int  blockId, iter, block[3];
	void initializeParticle(ptype x0, ptype y0, ptype z0, int tag);
	void track(ptype *W[5], UBspline_3d_d *splinep[5], ptype dt);
	void identifyBlock();
	void identifyBlockId();
	
	void updatePosition(ptype *W[5], UBspline_3d_d *splinep[5], ptype dt);	
	void collectInfo();
	void distributeInfo();
	void getParticleProp(UBspline_3d_d *splinep[5]);

private:
	
};

#endif
