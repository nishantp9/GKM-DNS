/***********************************************************************
 *                        6 MARCH 2014            
 *                      NISHANT PARASHAR         
 *                   DNS 3D Gas Kinetic Scheme    
 * --------------------------------------------------------------------*
 *                        MAIN ROUTINE
 * *********************************************************************/
 
#include <iostream>
#include <stdio.h>
#include <fstream>
#include "param.h"
#include "FlowField.h"

#ifdef UseCPU
	#include "CPUupdate.h"
#endif

#ifdef UseGPU
	#include "GPUupdate.h"
#endif

using namespace std;

void Iterate(Macro *macro);
void PrintW(Macro *macro, int i);
void PrintStats(Macro *macro);


/***********************************************************************
 * Main Routine
 ***********************************************************************/

int main(int argc, char* argv[])
{	
	Macro *macro = new Macro();

/* ----------------------------------------------------------*
 *  Initialize Flow field with DHIT conditions
 * ----------------------------------------------------------*/
	macro->Initialize_DHIT();
	macro->PeriodicBC();
	
	cout<<endl<<"------------GKM DNS-DHIT-----------"<<endl;
	cout<<"N    =  " << nc <<"x"<< nc <<"x"<< nc <<endl;
	cout<<"M    =  " << macro->Mt << endl;
	cout<<"k0   =  " << macro->k0 << endl;
		
/* ---------------------------*
 *  Update Flow Field
 * ---------------------------*/
	Iterate(macro);
	
	return 0;
				
}

/***********************************************************************
 * End of Main routine
 ***********************************************************************/


/* --------------------------------------------------------------------*
 * Evolve Using CPU, if using CPU else use GPU
 * --------------------------------------------------------------------*/

void Iterate(Macro *macro)
{
	#ifdef UseGPU
		cout<<"------------Reference TKE-----------"<<endl;
		cout<<endl<<"Time   TKE"<<endl;
		PrintStats(macro);
		cout<<"_______________________________________________________"<<endl;

		cout<<endl<<"-------Tracking Turbulence Statistics-------------"<<endl;
		cout<<endl<<"Time  	      TKE   		PRMS"<<endl;
		cout<<"________________________________________________________"<<endl;

		macro->dt = .0005;
		evolve(macro);
	#endif
	
/**************************************************************************/

	#ifdef UseCPU
		CPUevolve pu;
	
		cout<<"------------Reference TKE-------------"<<endl;
		cout<<endl<<"Time   TKE"<<endl;
		PrintStats(macro);
		cout<<"______________________________________"<<endl;
		// macro->timestep();
		macro->dt = .0005;
	
		cout<<endl<<"-------Tracking TKE decay rate--------"<<endl;
		cout<<endl<<"Time  		    TKE"<<endl;
		cout<<"_______________________________________"<<endl;

		while (macro->t < tmax)
		{
			pu.evolve(macro);
			//	macro->timestep();
			macro->t += macro->dt;				
			macro->PeriodicBC();
/**         For printing only TKE	**/
			PrintStats(macro);		

/**         For printing TKE, PRMS, TRMS, VRMS 		**/ 			
			//	pu.PrintStats(macro);		
			
		}
	
	#endif
}

/* --------------------------------------------------------------------*
 * Print W[i][:]
 * --------------------------------------------------------------------*/
void PrintW (Macro *macro, int l)
{
	FOR(j, nt)
	{	
		FOR(i, nt)
		{	
			FOR(k, nt)
				cout<<macro->W[l][I(j, i, k)]<<" ";
			
			cout<<endl;
		}
		cout<<endl;
	}
}

void PrintStats(Macro *macro)
{
	std::cout.precision(8);
	ptype TKE = 0;

	For(i, nt) 
		For(j, nt)
			For(k, nt)
			{	
				TKE = TKE + (pow(macro->W[1][I(i, j, k)],2) +
				             pow(macro->W[2][I(i, j, k)],2) +
				             pow(macro->W[3][I(i, j, k)],2)) /
				             (macro->W[0][I(i, j, k)]*pow(macro->c0*macro->Mt,2));	
			}

	TKE = TKE/Nc;

	std::cout<<macro->t/t0<<"       "<<TKE<<std::endl;		

}
/***********************************************************************
 ***********************************************************************/
