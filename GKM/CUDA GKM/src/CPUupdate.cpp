#include "CPUupdate.h"
#include "GKMflux.h"
#include <iostream>
using namespace std;

CPUevolve::CPUevolve()
{
	MAKE5(Fx);
	MAKE5(Fy);
	MAKE5(Fz);	
}

CPUevolve::~CPUevolve()
{
	KILL5(Fx);
	KILL5(Fy);
	KILL5(Fz);	
}


void CPUevolve::evolve(Macro *macro)
{
	derivsX(macro);
	derivsY(macro);
	derivsZ(macro);
	update(macro);
}

void CPUevolve::derivsX(Macro *macro)
{
	For(j, nt)
		FoR(i, nt)
		{	For(k, nt)
			{
				FOR(q, 5)
				{	
					WL[q] = macro->W[q][I(j, i, k)];
					WR[q] = macro->W[q][I(j, i+1, k)];
					W0[q] = .5*(WR[q] + WL[q]);
					
					WLN[q] = macro->W[q][I(j+1, i, k)];
					WRN[q] = macro->W[q][I(j+1, i+1, k)];
					W0N[q] = .5*(WRN[q] + WLN[q]);
										
					WLS[q] = macro->W[q][I(j-1, i, k)];
					WRS[q] = macro->W[q][I(j-1, i+1, k)];
					W0S[q] = .5*(WRS[q] + WLS[q]);
															
					WLF[q] = macro->W[q][I(j, i, k+1)];
					WRF[q] = macro->W[q][I(j, i+1, k+1)];
					W0F[q] = .5*(WRF[q] + WLF[q]);
										
					WLB[q] = macro->W[q][I(j, i, k-1)];
					WRB[q] = macro->W[q][I(j, i+1, k-1)];
					W0B[q] = .5*(WRB[q] + WLB[q]);
															
					DWx[q] = (WR[q]  - WL[q])/macro->dx;	
					DWy[q] = 0.5*(W0N[q] - W0S[q])/macro->dx;
					DWz[q] = 0.5*(W0F[q] - W0B[q])/macro->dx;

				}
				c2p();
//				macro->mu = macro->mu0*pow((Ti/macro->T0),1.5)*(macro->T0+110.4)/(Ti+110.4);

				macro->mu = macro->mu0*pow((Ti/macro->T0),.76);

			   flux(W0, DWx, DWy, DWz, macro->mu, macro->dt, macro->dx, F);
				
				FOR(r, 5)
					Fx[r][I(j, i, k)] = F[r];
			}

		}
}

void CPUevolve::derivsY(Macro *macro)
{
	ptype temp;
	FoR(j, nt)
		For(i, nt)
			For(k, nt)
			{
				FOR(q, 5)
				{	
					WL[q] = macro->W[q][I(j, i, k)];
					WR[q] = macro->W[q][I(j+1, i, k)];
      			    W0[q] = .5*(WR[q] + WL[q]);
					
					WLN[q] = macro->W[q][I(j, i-1, k)];
					WRN[q] = macro->W[q][I(j+1, i-1, k)];
					W0N[q] = .5*(WRN[q] + WLN[q]);
										
					WLS[q] = macro->W[q][I(j, i+1, k)];
					WRS[q] = macro->W[q][I(j+1, i+1, k)];
					W0S[q] = .5*(WRS[q] + WLS[q]);
															
					WLF[q] = macro->W[q][I(j, i, k+1)];
					WRF[q] = macro->W[q][I(j+1, i, k+1)];
					W0F[q] = .5*(WRF[q] + WLF[q]);
										
					WLB[q] = macro->W[q][I(j, i, k-1)];
					WRB[q] = macro->W[q][I(j+1, i, k-1)];
					W0B[q] = .5*(WRB[q] + WLB[q]);
															
					DWx[q] = (WR[q]  - WL[q])/macro->dx;	
					DWy[q] = 0.5*(W0S[q] - W0N[q])/macro->dx;
					DWz[q] = 0.5*(W0F[q] - W0B[q])/macro->dx;
					
				}
				c2p();
//				macro->mu = macro->mu0*pow((Ti/macro->T0),1.5)*(macro->T0+110.4)/(Ti+110.4);
				
				macro->mu = macro->mu0*pow((Ti/macro->T0),.76);

				temp = W0[2];   W0[2] = W0[1];   W0[1] =  temp;

				temp = DWx[2]; DWx[2] = DWx[1]; DWx[1] =  temp;
				temp = DWy[2]; DWy[2] = DWy[1]; DWy[1] =  temp;
				temp = DWz[2]; DWz[2] = DWz[1]; DWz[1] =  temp;


			   flux(W0, DWx, DWy, DWz, macro->mu, macro->dt, macro->dx, F);
				
				FOR(r, 5)
					Fy[r][I(j, i, k)] = F[r];
				
				temp  =  Fy[2][I(j, i, k)]; 
				Fy[2][I(j, i, k)] =  Fy[1][I(j, i, k)];
				Fy[1][I(j, i, k)] =  temp;
			}

}	

void CPUevolve::derivsZ(Macro *macro)
{
	ptype temp;
	For(j, nt)
		For(i, nt)
		{	FoR(k, nt)
			{
				FOR(q, 5)
				{	
					WL[q] = macro->W[q][I(j, i, k)];
					WR[q] = macro->W[q][I(j, i, k+1)];
      			    W0[q] = .5*(WR[q] + WL[q]);
					
					WLN[q] = macro->W[q][I(j, i+1, k)];
					WRN[q] = macro->W[q][I(j, i+1, k+1)];
					W0N[q] = .5*(WRN[q] + WLN[q]);
										
					WLS[q] = macro->W[q][I(j, i-1, k)];
					WRS[q] = macro->W[q][I(j, i-1, k+1)];
					W0S[q] = .5*(WRS[q] + WLS[q]);
															
					WLF[q] = macro->W[q][I(j+1, i, k)];
					WRF[q] = macro->W[q][I(j+1, i, k+1)];
					W0F[q] = .5*(WRF[q] + WLF[q]);
										
					WLB[q] = macro->W[q][I(j-1, i, k)];
					WRB[q] = macro->W[q][I(j-1, i, k+1)];
					W0B[q] = .5*(WRB[q] + WLB[q]);
															
					DWx[q] = (WR[q]  - WL[q])/macro->dx;	
					DWy[q] = 0.5*(W0N[q] - W0S[q])/macro->dx;
					DWz[q] = 0.5*(W0F[q] - W0B[q])/macro->dx;

				}
				c2p();
//				macro->mu = macro->mu0*pow((Ti/macro->T0),1.5)*(macro->T0+110.4)/(Ti+110.4);

				macro->mu = macro->mu0*pow((Ti/macro->T0),.76);

				temp = W0[3]; W0[3]  = W0[1];  W0[1]  =  temp;

				temp = DWx[3]; DWx[3] = DWx[1]; DWx[1] =  temp;
				temp = DWy[3]; DWy[3] = DWy[1]; DWy[1] =  temp;
				temp = DWz[3]; DWz[3] = DWz[1]; DWz[1] =  temp;
				
				
				flux(W0, DWx, DWy, DWz, macro->mu, macro->dt, macro->dx, F);				
				
				FOR(r, 5)
					Fz[r][I(j, i, k)] = F[r];

				temp  =  Fz[3][I(j, i, k)]; 
				Fz[3][I(j, i, k)] =  Fz[1][I(j, i, k)];
				Fz[1][I(j, i, k)] =  temp;

			}
		}
	
}

void CPUevolve::update(Macro *macro)
{   	
	For(j, nt)
		For(i, nt)
			For(k, nt)
				FOR(q, 5)
				{					
					macro->W[q][I(j, i, k)]  =  
					macro->W[q][I(j, i, k)]  -
					(1/macro->dx) * (Fx[q][I(j, i, k)]  - Fx[q][I(j, i-1, k)] 
										 +
									     Fy[q][I(j, i, k)]  - Fy[q][I(j-1, i, k)] 
									     +
									     Fz[q][I(j, i, k)]  - Fz[q][I(j, i, k-1)]);
            }
           
}

void CPUevolve::PrintStats(Macro *macro)
{
	std::cout.precision(8);
	ptype TKE = 0;
	ptype Pmean = 0;	ptype Prms = 0;
	ptype Tmean = 0;	ptype Trms = 0;
	ptype Vmean = 0;	ptype Vrms = 0;
	
	ptype umean = 0;	ptype urms = 0;
	ptype vmean = 0;	ptype vrms = 0;
	ptype wmean = 0;	ptype wrms = 0;
	
	

	For(i, nt) 
		For(j, nt)
			For(k, nt)
				{
					FOR(q, 5)
						W0[q] = macro->W[q][I(i, j, k)]; 					
					
					c2p();
						
//					TKE = TKE + ( Vsqr / pow(macro->c0*macro->Mt,2));
								
					Pmean = Pmean + Pi;
					Tmean = Tmean + Ti;
					Vmean = Vmean + 1/deni; 				

					umean = umean + ui; 				
					vmean = vmean + vi; 				
					wmean = wmean + wi; 				

				}

	Pmean /= Nc;
	Tmean /= Nc;
	Vmean /= Nc;
	
	umean /= Nc;
	vmean /= Nc;
	wmean /= Nc;


	For(i, nt) 
		For(j, nt)
			For(k, nt)
				{

					FOR(q, 5)
						W0[q] = macro->W[q][I(i, j, k)];
					
					c2p();
										   
				   Prms = Prms + pow((Pmean-Pi),2);
				   Trms = Trms + pow((Tmean-Ti),2);
				   Vrms = Vrms + pow((Vmean-1/deni),2);				

				   urms = urms + pow((umean-ui),2);				
				   vrms = vrms + pow((vmean-vi),2);				
				   wrms = wrms + pow((wmean-wi),2);				
				}

	Prms = sqrt(Prms/Nc) / (gam*macro->p0*macro->Mt*macro->Mt);
	Trms = sqrt(Trms/Nc) / ((gam-1)*macro->T0*macro->Mt*macro->Mt);
	Vrms = sqrt(Vrms/Nc) / (macro->Mt*macro->Mt / macro->den0);

	urms = sqrt(urms/Nc);
	vrms = sqrt(vrms/Nc);
	wrms = sqrt(wrms/Nc);

	TKE = ( urms*urms + vrms*vrms + wrms*wrms ) / pow(macro->c0*macro->Mt,2);
	
	std::cout<<macro->t/t0<<"  "<<TKE<<"  "<<Prms<<"  "<<Trms<<"  "<<Vrms<<std::endl;		
}	



void CPUevolve::c2p()
{
	deni = W0[0];
	ui   = W0[1]/deni;
	vi   = W0[2]/deni;
	wi   = W0[3]/deni;
	Vsqr = (ui*ui + vi*vi + wi*wi);
	Pi   = (deni*(gam-1))*(W0[4]/deni - .5*Vsqr);
	Ti   = Pi/deni;
}
 
void CPUevolve::PrintF (ptype *F[], int l)
{

	For(j, nt)
	{	
		For(i, nt)
		{	
			For(k, nt)
				cout<<F[l][I(j, i, k)]<<" ";
			
			cout<<endl;
		}
		cout<<endl;
	}
}
