/***************************************************************
 * FlowField Definitions
 ***************************************************************/

#include "FlowField.h"
#include <fstream>
#include <iostream>
#include <stdio.h>
using namespace std;
/* ------------------------------------------------------*
 *  Overridden Constructor
 * ------------------------------------------------------*/
Macro::Macro()
{
	ifstream input(DHIT_Binary, ios::in | ios::binary);
	int n;
	input.read((char*) &n, sizeof(int));
	input.close();
	
	::nc = n;
	::nt = nc+2;
	::Nc = nc*nc*nc;
	::Nt = nt*nt*nt;
	t = 0;
	MAKE5(W);

	MAKE_(u);
	MAKE_(v);
	MAKE_(w);
	MAKE_(c);	
	MAKE_(den);	
	MAKE_(p);
	MAKE_(M);		
}

/* -------------------------------*
 *  Overridden Destructor
 * -------------------------------*/
Macro::~Macro()
{
	KILL5(W);

	KILL_(u);
	KILL_(v);
	KILL_(w);
	KILL_(c);	
	KILL_(den);	
	KILL_(p);	
	KILL_(M);	
}

/* ------------------------------------------------------*
 * Conserative --> Primitive 3D
 * ------------------------------------------------------*/ 
void Macro::c2p3D()
{
	int I, J;
	FOR(j, nc)
		FOR(i, nc)
			FOR(k, nc)
			{
				I = Ii(j, i, k);
				J = Ic(j, i, k);
				den[J] = W[0][I];
				u[J] = W[1][I]/W[0][I];
				v[J] = W[2][I]/W[0][I];
				w[J] = W[3][I]/W[0][I];
				Vsqr = (u[J]*u[J] + v[J]*v[J] + w[J]*w[J]);
			
				p[J] = W[4][I]*(gam-1) - 0.5*Vsqr*W[0][I]*(gam-1);
				c[J] = sqrt( gam*p[J] / W[0][I] );
			}
}	

/* -----------------------------------------------------*
 * Primitive to Conservative 3D
 * -----------------------------------------------------*/ 	
void Macro::p2c3D()
{	
	int I, J;
	FOR(j, nc)
		FOR(i, nc)
			FOR(k, nc)
			{	
				I = Ii(j, i, k);
				J = Ic(j, i, k);
				
				Vsqr = (u[J]*u[J] + v[J]*v[J] + w[J]*w[J]);
				W[0][I] = den[J];
				W[1][I] = u[J]*den[J];
				W[2][I] = v[J]*den[J];
				W[3][I] = w[J]*den[J];
				W[4][I] = den[J]*(p[J]/(den[J]*(gam-1.0)) + .5*Vsqr);	   
			}
}

/* -----------------------------------------------------*
 * CFL timestep calculation
 * -----------------------------------------------------*/  
void Macro::timestep()
{
/**	
	c2p3D();
	maxU = max(u);
	maxC = max(c);
	ptype invRe = (den0*maxU*dx) / mu;
	dt = cfl*dx/((maxU + maxC)*(1 + .5*invRe));	
**/

	ptype den, Ux, Uy, Uz, U, P, valC;
	maxU = 0;
	maxC = 0;
	
	For(j, nt)
		For(i, nt)
			For(k, nt)
			{ 
				den = W[0][I(j, i, k)];
				Ux  = W[1][I(j, i, k)] / den;
				Uy  = W[2][I(j, i, k)] / den;
				Uz  = W[3][I(j, i, k)] / den;
				
				U   = sqrt(Ux*Ux+Uy*Uy+Uz*Uz);
				
				P = (den*(gam-1))*(W[4][I(j, i, k)] / den -
								.5*Ux*Ux-.5*Uy*Uy-.5*Uz*Uz);
		
				if(U > maxU) {
					maxU = U;
				}
				
				valC = sqrt(gam*P/den);
				if(valC > maxC) {
					maxC = valC;
				}

			}

	ptype invRe = 2 * mu / (den0*maxU*dx);
	dt = cfl*dx/((maxU + maxC)*(1 + invRe));	


/**	ptype den, Ux, Uy, Uz, U, P, C, invRe, dt_min = dx;
	
	For(j, nt)
		For(i, nt)
			For(k, nt)
			{ 
				den = W[0][I(j, i, k)];
				Ux  = W[1][I(j, i, k)] / den;
				Uy  = W[2][I(j, i, k)] / den;
				Uz  = W[3][I(j, i, k)] / den;
				
				U   = sqrt(Ux*Ux+Uy*Uy+Uz*Uz);
				
				
				P = (den*(gam-1))*(W[4][I(j, i, k)] / den -
								.5*Ux*Ux-.5*Uy*Uy-.5*Uz*Uz);
		
				C = sqrt(gam*P/den);
				
				invRe = 2 * mu / (den*U*dx);
				
				dt = cfl*dx/((U + C)*(1 + invRe));
				
				if(dt < dt_min)
				{
					dt_min = dt;
				}		
			}
**/
}

/*------------------------------------------------------*
 * Initialization with DHIT initial conditions
 * -----------------------------------------------------*/

void Macro::Initialize_DHIT()
{
	int n;
	ifstream input(DHIT_Binary, ios::in | ios::binary);
	input.read((char*) &n, sizeof(int));
	
	size_t size = sizeof(ptype);
	input.read((char*) &den0, size);
	input.read((char*) &Mt, size);
	input.read((char*) &Re, size);
	input.read((char*) &k0, size);
	input.read((char*) &urms0, size);
	input.read((char*) &mu0, size);
	input.read((char*) &c0, size);
	input.read((char*) &p0, size);
	ptype t;
	input.read((char*) &t, size);
	::t0 = t;
	dx  =  2.0*pi/nc;
   T0  =  p0/den0;

/* -------------------------------------------------------------------------------------*
 * storing primitive initial conds to conservative variable
 * -------------------------------------------------------------------------------------*/
	ptype *u0, *v0, *w0;
	u0 = new ptype[Nc];
	v0 = new ptype[Nc];
	w0 = new ptype[Nc];
	
	input.read((char*) u0, size*Nc);
	input.read((char*) v0, size*Nc);
	input.read((char*) w0, size*Nc);
	
	FOR(j, nc)
		FOR(i, nc)
			FOR(k, nc)
			{
				// we adopt the following index(j,i,k) style for linearization
				int I = Ii(j, i, k);
				
				// Matlab style of 3-D array linearization --> our style equiv. index(k,i,j) 
				int J = k*nc*nc + i*nc + j; 	
				
				Vsqr = (u0[J]*u0[J] + v0[J]*v0[J] + w0[J]*w0[J]);
				W[0][I] = den0;
				W[1][I] = u0[J]*den0;
				W[2][I] = v0[J]*den0;
				W[3][I] = w0[J]*den0;
				W[4][I] = den0*(p0/(den0*(gam-1.0)) + .5*Vsqr);
			}
			
	delete [] u0;
	delete [] v0;
	delete [] w0;
	
	input.close();
	
}

/* ------------------------------------------------------*
 * Periodic boundary conditions
 * ------------------------------------------------------*/
void Macro::PeriodicBC()
{
/** X-direction Periodic BC **/ 
	FOR(i, 5)
		FOR(j, nt)
			FOR(k, nt)
			{
				W[i][I(j, 0, k)]     = 	W[i][I(j, nt-2, k)];
				W[i][I(j, nt-1, k)]  = 	W[i][I(j, 1, k)];
			}

/** Y-direction Periodic BC **/
	FOR(i, 5)
		FOR(j, nt)
			FOR(k, nt)
			{
				W[i][I(0, j, k)]     = 	W[i][I(nt-2, j, k)];
				W[i][I(nt-1, j, k)]  = 	W[i][I(1, j, k)];
			}

/** Z-direction Periodic BC **/
	FOR(i, 5)
		FOR(j, nt)
			FOR(k, nt)
			{
				W[i][I(j, k, 0)]     = 	W[i][I(j, k, nt-2)];
				W[i][I(j, k, nt-1)]  = 	W[i][I(j, k, 1)];
			}	

}

void Macro::SaveVelocityField()
{
	c2p3D();
/**	
	fstream u_writer("InitialConditions/ u.bin", std::ios::out | std::ios::binary);
	fstream v_writer("InitialConditions/ v.bin", std::ios::out | std::ios::binary);
	fstream w_writer("InitialConditions/ w.bin", std::ios::out | std::ios::binary);
	
	size_t size = sizeof(ptype) * Nc; 
	u_writer.write((char*) u, size);
	v_writer.write((char*) v, size);
	w_writer.write((char*) w, size);
	
	u_writer.flush();
	v_writer.flush();
	w_writer.flush();
	
	u_writer.close();
	v_writer.close();
	w_writer.close();
**/

	std::ofstream outfile1 ("InitialConditions/u.dat");
	std::ofstream outfile2 ("InitialConditions/v.dat");
	std::ofstream outfile3 ("InitialConditions/w.dat");
	std::ofstream outfile4 ("InitialConditions/M.dat");
	
	FOR(i, Nc)
	{
		M[i] = sqrt(u[i]*u[i] + v[i]*v[i] + w[i]*w[i]) / c[i];
		outfile1 <<u[i]<<endl;
		outfile2 <<v[i]<<endl;
		outfile3 <<w[i]<<endl;
		outfile4 <<M[i]<<endl;
	}
	outfile1.close();
	outfile2.close();
	outfile3.close();
	outfile4.close();
}

void Macro::SaveVelocityField2()
{
	c2p3D();
/**	
	fstream u_writer("InitialConditions/ u.bin", std::ios::out | std::ios::binary);
	fstream v_writer("InitialConditions/ v.bin", std::ios::out | std::ios::binary);
	fstream w_writer("InitialConditions/ w.bin", std::ios::out | std::ios::binary);
	
	size_t size = sizeof(ptype) * Nc; 
	u_writer.write((char*) u, size);
	v_writer.write((char*) v, size);
	w_writer.write((char*) w, size);
	
	u_writer.flush();
	v_writer.flush();
	w_writer.flush();
	
	u_writer.close();
	v_writer.close();
	w_writer.close();
**/
	ptype T;
	FILE * pFile;
	pFile = fopen ("DNS.dat","w");


	FOR(i, Nc)
	{
		T = p[i]/(287*den[i]);
		
		fprintf(pFile, "%.15f  %.15f  %.15f  %.15f  %.15f  %.15f\n", u[i], v[i], w[i], den[i], p[i], T );
		
	}
	
	fclose(pFile);
}

void Macro::SaveVelocityField3()
{
	c2p3D();
/**	
	fstream u_writer("InitialConditions/ u.bin", std::ios::out | std::ios::binary);
	fstream v_writer("InitialConditions/ v.bin", std::ios::out | std::ios::binary);
	fstream w_writer("InitialConditions/ w.bin", std::ios::out | std::ios::binary);
	
	size_t size = sizeof(ptype) * Nc; 
	u_writer.write((char*) u, size);
	v_writer.write((char*) v, size);
	w_writer.write((char*) w, size);
	
	u_writer.flush();
	v_writer.flush();
	w_writer.flush();
	
	u_writer.close();
	v_writer.close();
	w_writer.close();
**/
	ptype T;
	FILE * pFile;
	pFile = fopen ("DNS2.dat","w");


	FOR(i, Nc)
	{
		T = p[i]/(287*den[i]);
		
		fprintf(pFile, "%.15f  %.15f  %.15f  %.15f  %.15f  %.15f\n", u[i], v[i], w[i], den[i], p[i], T );
		
	}
	
	fclose(pFile);
}
/***************************************************************
 * --------------------------END-------------------------------*/
