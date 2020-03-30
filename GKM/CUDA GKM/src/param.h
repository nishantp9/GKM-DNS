/****************************************************************************************
 * Stores the basic problem parameters
 * -------------------------------------------------------------------------------------*
 ****************************************************************************************/

#ifndef PARAMHEADERDEF
#define PARAMHEADERDEF

#include <cmath>

extern int nc, nt, Nc, Nt;

#define UseGPU
//#define UseCPU
#define	pi 	atan(1)*4
#define  eps   pow(10,-10)
#define  K  	2
#define  gam   7.0/5.0
#define  PRN   1
#define	cfl   0.1

#define DHIT_Binary "InitialConditions/Init.dat"

//typedef double ptype;	//precision type
typedef float double;	//precision type

/* ------------------------------------------*
 * short-hand macros
 * ------------------------------------------*/
#define FOR(i, n) for(int i=0; i<n; i++)
#define For(i, n) for(int i=1; i<n-1; i++)
#define FoR(i, n) for(int i=0; i<n-1; i++)

// we adopt the following index(i,j,k) style for linearization of 3-D array
#define I(i, j, k)  (i)*nt*nt + (j)*nt + (k) 
// note :: brackets around i, j & k are put on purpose
// so that it can handle args with arithmetic operands e.g I(nt-4, j+1 ,k)	

// Index def used for reading Initial conditions (nc x nc x nc) 
#define Ii(i, j, k) (i+1)*nt*nt + (j+1)*nt + (k+1)
#define Ic(i, j, k) (i)*nc*nc + (j)*nc + (k)

/* -------------------------------------------------------------------------------------*
 * macros for declaring and deleting dynamically allocated memory
 * -------------------------------------------------------------------------------------*/
#define MAKE(var)  var = new ptype[Nt]
#define KILL(var)  delete [] var
#define MAKE5(var) MAKE(var[0]); MAKE(var[1]); MAKE(var[2]); MAKE(var[3]); MAKE(var[4]);
#define KILL5(var) KILL(var[0]); KILL(var[1]); KILL(var[2]); KILL(var[3]); KILL(var[4]); 

// Save Files (without ghost nodes) for Post-Processing in  Matlab
#define MAKE_(var)  var = new ptype[Nc]
#define KILL_(var)  delete [] var

extern ptype t0;  
#define  tmax  3*t0

// 1) max template
template <typename T> T max(T *A)
{
	T Max = A[0];
	FOR(i, Nt)
	{			
		if(A[i] > Max)
			Max = A[i];
	}
	return Max;
}

// 2) sign template
template <typename T> int sign(T s)
{
	int sign;
	if(s>0)
	sign=1;
	else if(s<0)
	sign=-1;
	else if(s==0)
	sign=0;
	return sign;
}

#endif

/****************************************************************************************
 * -------------------------------------------------------------------------------------*/
