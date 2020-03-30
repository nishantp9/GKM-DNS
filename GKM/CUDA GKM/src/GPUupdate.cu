/**************************************************************************************************************************
 * Updating Flow-Field on GPU --------------------------------------------------------------------------------------------*
 * -----------------------------------------------------------------------------------------------------------------------*
 **************************************************************************************************************************/

#include "GPUupdate.h"
#include <iostream>

using namespace std;

#define Bsize 8
#define Bsz3D 512 

#define CUMAKE(var)  cudaMalloc((void**) &var, sizeof(ptype)*Nt)
#define CUKILL(var)  cudaFree(var)

#define CUMAKE5(var) CUMAKE(var[0]); CUMAKE(var[1]); CUMAKE(var[2]); CUMAKE(var[3]); CUMAKE(var[4]);
#define CUKILL5(var) CUKILL(var[0]); CUKILL(var[1]); CUKILL(var[2]); CUKILL(var[3]); CUKILL(var[4]);

#define HtoD(dv,v) cudaMemcpy(dv, v, sizeof(ptype)*Nt, cudaMemcpyHostToDevice)
#define DtoH(dv,v) cudaMemcpy(v, dv, sizeof(ptype)*Nt, cudaMemcpyDeviceToHost)

#define HtoD5(dv,v) HtoD(dv[0],v[0]); HtoD(dv[1],v[1]); HtoD(dv[2],v[2]); HtoD(dv[3],v[3]); HtoD(dv[4],v[4]); 
#define DtoH5(dv,v) DtoH(dv[0],v[0]); DtoH(dv[1],v[1]); DtoH(dv[2],v[2]); DtoH(dv[3],v[3]); DtoH(dv[4],v[4]);
 
#define CUMAKE_(var)  cudaMalloc((void**) &var, sizeof(ptype)*Nc)
#define CUKILL_(var)  cudaFree(var)

#define CUMAKE5_(var) CUMAKE_(var[0]); CUMAKE_(var[1]); CUMAKE_(var[2]); CUMAKE_(var[3]); CUMAKE_(var[4]);
#define CUKILL5_(var) CUKILL_(var[0]); CUKILL_(var[1]); CUKILL_(var[2]); CUKILL_(var[3]); CUKILL_(var[4]);

#define HtoD_(dv,v) cudaMemcpy(dv, v, sizeof(ptype)*Nc, cudaMemcpyHostToDevice)
#define DtoH_(dv,v) cudaMemcpy(v, dv, sizeof(ptype)*Nc, cudaMemcpyDeviceToHost)

#define HtoD5_(dv,v) HtoD_(dv[0],v[0]); HtoD_(dv[1],v[1]); HtoD_(dv[2],v[2]); HtoD_(dv[3],v[3]); HtoD_(dv[4],v[4]); 
#define DtoH5_(dv,v) DtoH_(dv[0],v[0]); DtoH_(dv[1],v[1]); DtoH_(dv[2],v[2]); DtoH_(dv[3],v[3]); DtoH_(dv[4],v[4]);



/**************************************************************************************************************************/


__constant__ int d_nt, d_nc, d_Nt, d_Nc;

__constant__ ptype dx, dt, mu0, T0, GAM;

__global__ void derivsX_Kernel(ptype *d_W, ptype *d_We, ptype *d_DWx, ptype *d_DWy, ptype *d_DWz);

__global__ void derivsY_Kernel(ptype *d_W, ptype *d_We, ptype *d_DWx, ptype *d_DWy, ptype *d_DWz);

__global__ void derivsZ_Kernel(ptype *d_W, ptype *d_We, ptype *d_DWx, ptype *d_DWy, ptype *d_DWz);

__global__ void flip_Kernel(ptype *d_W1, ptype *d_W2, ptype *d_DW1x, ptype *d_DW2x, ptype *d_DW1y,
                                                      ptype *d_DW2y, ptype *d_DW1z, ptype *d_DW2z);
																		 
__global__ void flipBack_Kernel(ptype *d_F1, ptype *d_F2);

__global__ void flux(ptype *d_W0, ptype *d_W1, ptype *d_W2, ptype *d_W3, ptype *d_W4,
				    ptype *d_DW0x, ptype *d_DW1x, ptype *d_DW2x, ptype *d_DW3x, ptype *d_DW4x,
				    ptype *d_DW0y, ptype *d_DW1y, ptype *d_DW2y, ptype *d_DW3y, ptype *d_DW4y,
			       ptype *d_DW0z, ptype *d_DW1z, ptype *d_DW2z, ptype *d_DW3z, ptype *d_DW4z,
				    ptype *d_F0, ptype *d_F1, ptype *d_F2, ptype *d_F3, ptype *d_F4, ptype *d_tau);


/*__global__ void update_Kernel (ptype *d_W0, ptype *d_W1, ptype *d_W2, ptype *d_W3, ptype *d_W4,
							   ptype *d_F0x, ptype *d_F1x, ptype *d_F2x, ptype *d_F3x, ptype *d_F4x,
						      ptype *d_F0y, ptype *d_F1y, ptype *d_F2y, ptype *d_F3y, ptype *d_F4y,
							   ptype *d_F0z, ptype *d_F1z, ptype *d_F2z, ptype *d_F3z, ptype *d_F4z);								  
*/

__global__ void update_Kernel (ptype *d_W, ptype *d_Fx, ptype *d_Fy, ptype *d_Fz);								  



__global__ void PeriodicBC_Kernel (ptype *d_W0, ptype *d_W1, ptype *d_W2, ptype *d_W3, ptype *d_W4);

__global__ void c2p3D_kernel (ptype *W0, ptype *W1, ptype *W2, ptype *W3, ptype *W4,
										ptype *u,  ptype *v,  ptype *w,  ptype *p,  ptype *r, ptype *k);
										 
__global__ void KE_kernel (ptype *KE, ptype *KE_reduced);

__global__ void PfluctuationSqr_Kernel(ptype *d_PflucSqr, ptype *d_p, ptype *d_Pmean);
										  
__device__ void d_slopesolver(ptype b[5], ptype U[3], ptype lam, ptype a[5]);

__device__ void d_mult(ptype M[5*2], int k, ptype C, ptype B[5]);

__device__ void d_Mult(ptype M[5*5], ptype C, ptype A[5], ptype B[5]);

__device__ void d_add(ptype A[5], ptype B[5]);

__device__ void d_MCal(ptype M[5],ptype If[3*7],ptype Ie2,ptype Ie4, int k,int l,int m, ptype ax[5]);
/***************************************************************************************************************************
 * ------------------------------------------------------------------------------------------------------------------------*
 ***************************************************************************************************************************/

void evolve(Macro* macro)
{
	cudaSetDevice(1);

	ptype *d_W0,   *d_W1,   *d_W2,   *d_W3,   *d_W4;
	ptype *d_W0e,  *d_W1e,  *d_W2e,  *d_W3e,  *d_W4e; 
	
	ptype *d_F0x,  *d_F1x,  *d_F2x,  *d_F3x,  *d_F4x;
	ptype *d_F0y,  *d_F1y,  *d_F2y,  *d_F3y,  *d_F4y;
	ptype *d_F0z,  *d_F1z,  *d_F2z,  *d_F3z,  *d_F4z;
		
	ptype *d_DW0x, *d_DW1x, *d_DW2x, *d_DW3x, *d_DW4x;
	ptype *d_DW0y, *d_DW1y, *d_DW2y, *d_DW3y, *d_DW4y;
	ptype *d_DW0z, *d_DW1z, *d_DW2z, *d_DW3z, *d_DW4z;
	
	ptype *d_tau;
		
	CUMAKE(d_W0);	CUMAKE(d_W1);	CUMAKE(d_W2);	CUMAKE(d_W3);	CUMAKE(d_W4);
	CUMAKE(d_W0e);	CUMAKE(d_W1e);	CUMAKE(d_W2e);	CUMAKE(d_W3e);	CUMAKE(d_W4e);
		
	CUMAKE(d_F0x);	CUMAKE(d_F1x);	CUMAKE(d_F2x);	CUMAKE(d_F3x);	CUMAKE(d_F4x);
	CUMAKE(d_F0y);	CUMAKE(d_F1y);	CUMAKE(d_F2y);	CUMAKE(d_F3y);	CUMAKE(d_F4y);
	CUMAKE(d_F0z);	CUMAKE(d_F1z);	CUMAKE(d_F2z);	CUMAKE(d_F3z);	CUMAKE(d_F4z);
		
	CUMAKE(d_DW0x);	CUMAKE(d_DW1x);	CUMAKE(d_DW2x);	CUMAKE(d_DW3x);	CUMAKE(d_DW4x);
	CUMAKE(d_DW0y);	CUMAKE(d_DW1y);	CUMAKE(d_DW2y);	CUMAKE(d_DW3y);	CUMAKE(d_DW4y);
	CUMAKE(d_DW0z);	CUMAKE(d_DW1z);	CUMAKE(d_DW2z);	CUMAKE(d_DW3z);	CUMAKE(d_DW4z);
	
	CUMAKE(d_tau);		
	
	ptype *d_r, *d_u, *d_v, *d_w, *d_p, *d_k, *d_k_reduced, *k_reduced, *d_PflucSqr, *d_Prms, *d_Pmean;
	
	CUMAKE_(d_r);	CUMAKE_(d_u);	CUMAKE_(d_v);	CUMAKE_(d_w);	CUMAKE_(d_p);	CUMAKE_(d_k);	
	CUMAKE_(d_PflucSqr);	
	cudaMalloc((void**) &d_Prms, sizeof(ptype));
	cudaMalloc((void**) &d_Pmean, sizeof(ptype));
	
	int sz_k = Nc / (Bsz3D<<1) + 1 ;

    	
	k_reduced = new ptype[sz_k];
	cudaMalloc((void**) &d_k_reduced, sizeof(ptype)*sz_k);
	
//	cudaMemset(d_F0z, 0, nt*sizeof(ptype));	
	
/** Host to Device Memcpy **/	
	HtoD(d_W0, macro->W[0]); 
	HtoD(d_W1, macro->W[1]); 
	HtoD(d_W2, macro->W[2]); 
	HtoD(d_W3, macro->W[3]); 
	HtoD(d_W4, macro->W[4]); 
		
	cudaThreadSynchronize();
			
	ptype Dx = macro->dx;
	ptype Dt = macro->dt;	
	ptype Mu = macro->mu0;
	ptype T	= macro->T0;
	ptype Ga = gam;
	
	cudaMemcpyToSymbol(d_nt, &nt, sizeof(int));
	cudaMemcpyToSymbol(d_nc, &nc, sizeof(int));
	cudaMemcpyToSymbol(d_Nt, &Nt, sizeof(int));
	cudaMemcpyToSymbol(d_Nc, &Nc, sizeof(int));
	cudaMemcpyToSymbol(dx  , &Dx, sizeof(ptype));
	cudaMemcpyToSymbol(dt  , &Dt, sizeof(ptype));
	cudaMemcpyToSymbol(mu0 , &Mu, sizeof(ptype));
	cudaMemcpyToSymbol(GAM , &Ga, sizeof(ptype));
	cudaMemcpyToSymbol(T0  ,  &T, sizeof(ptype));	

	ptype *F[5];

	MAKE(F[0]);	MAKE(F[1]);	MAKE(F[2]);	MAKE(F[3]);	MAKE(F[4]);	
	
	int m = 0;	
		
/***************************************************************************************************************************
 * 
 * ------------------------------------------------------------------------------------------------------------------------*/	
// Number of iterations
	int NumItr = tmax/macro->dt;
	FOR(h, NumItr) {
	
/*	HtoD(d_W0, macro->W[0]); 
	HtoD(d_W1, macro->W[1]); 
	HtoD(d_W2, macro->W[2]); 
	HtoD(d_W3, macro->W[3]); 
	HtoD(d_W4, macro->W[4]); 	
	cudaThreadSynchronize();
*/		
	int blocks = (nt-1)/(Bsize-2) + 1;
	dim3 dimBlock(Bsize, Bsize, Bsize);
	dim3 dimGrid (blocks, blocks, blocks);
		
	derivsX_Kernel<<<dimGrid, dimBlock>>>(d_W0, d_W0e, d_DW0x, d_DW0y, d_DW0z);
	derivsX_Kernel<<<dimGrid, dimBlock>>>(d_W1, d_W1e, d_DW1x, d_DW1y, d_DW1z);
	derivsX_Kernel<<<dimGrid, dimBlock>>>(d_W2, d_W2e, d_DW2x, d_DW2y, d_DW2z);
	derivsX_Kernel<<<dimGrid, dimBlock>>>(d_W3, d_W3e, d_DW3x, d_DW3y, d_DW3z);
	derivsX_Kernel<<<dimGrid, dimBlock>>>(d_W4, d_W4e, d_DW4x, d_DW4y, d_DW4z);
	cudaThreadSynchronize();
		

	flux<<<dimGrid, dimBlock>>>(d_W0e, d_W1e, d_W2e, d_W3e, d_W4e, d_DW0x, d_DW1x, d_DW2x, d_DW3x, d_DW4x, d_DW0y,
															d_DW1y,d_DW2y, d_DW3y, d_DW4y, d_DW0z, d_DW1z, d_DW2z,
															d_DW3z, d_DW4z, d_F0x, d_F1x, d_F2x, d_F3x, d_F4x, d_tau);


/**********************
 * CHECK with CPU code
 **********************/
	
/*	DtoH(d_F0x, F[0]);
	cudaThreadSynchronize();	
	DtoH(d_F1x, F[1]);
	cudaThreadSynchronize(); 	
	DtoH(d_F2x, F[2]);
	cudaThreadSynchronize();
	DtoH(d_F3x, F[3]);
	cudaThreadSynchronize();
	DtoH(d_F4x, F[4]);
	cudaThreadSynchronize();

	For(i, nt)
	{	
		For(j, nt)
		{	
			For(k, nt)
				cout<<F[0][i*nt*nt + j*nt + k]<<" ";
			
			cout<<endl;
		}
		cout<<endl;
	}
*/
/**************************************************************************************************************************
 * -----------------------------------------------------------------------------------------------------------------------*/      
	
	derivsY_Kernel<<<dimGrid, dimBlock>>>(d_W0, d_W0e, d_DW0x, d_DW0y, d_DW0z);
	derivsY_Kernel<<<dimGrid, dimBlock>>>(d_W1, d_W1e, d_DW1x, d_DW1y, d_DW1z);
	derivsY_Kernel<<<dimGrid, dimBlock>>>(d_W2, d_W2e, d_DW2x, d_DW2y, d_DW2z);
	derivsY_Kernel<<<dimGrid, dimBlock>>>(d_W3, d_W3e, d_DW3x, d_DW3y, d_DW3z);
	derivsY_Kernel<<<dimGrid, dimBlock>>>(d_W4, d_W4e, d_DW4x, d_DW4y, d_DW4z);
	cudaThreadSynchronize();

	flip_Kernel<<<dimGrid, dimBlock>>>(d_W1e, d_W2e, d_DW1x, d_DW2x, d_DW1y, d_DW2y, d_DW1z, d_DW2z);
	cudaThreadSynchronize();

	flux<<<dimGrid, dimBlock>>>(d_W0e, d_W1e, d_W2e, d_W3e, d_W4e, d_DW0x, d_DW1x, d_DW2x, d_DW3x, d_DW4x, d_DW0y,
															d_DW1y,d_DW2y, d_DW3y, d_DW4y, d_DW0z, d_DW1z, d_DW2z,
															d_DW3z, d_DW4z, d_F0y, d_F1y, d_F2y, d_F3y, d_F4y, d_tau);
   cudaThreadSynchronize();
/*
	DtoH(d_F0y, F[0]);
	cudaThreadSynchronize();	
	DtoH(d_F1y, F[1]);
	cudaThreadSynchronize(); 	
	DtoH(d_F2y, F[2]);
	cudaThreadSynchronize();
	DtoH(d_F3y, F[3]);
	cudaThreadSynchronize();
	DtoH(d_F4y, F[4]);
	cudaThreadSynchronize();

	FOR(i, nt)
	{	
		FOR(j, nt)
		{	
			FOR(k, nt)
				cout<<F[3][i*nt*nt + j*nt + k]<<" ";
			
			cout<<endl;
		}
		cout<<endl;
	}
*/

   flipBack_Kernel<<<dimGrid, dimBlock>>>(d_F1y, d_F2y);
   cudaThreadSynchronize();


/**********************
 * CHECK with CPU code
 **********************/
/*	DtoH(d_F0y, F[0]);
	cudaThreadSynchronize();	
	DtoH(d_F1y, F[1]); 
	cudaThreadSynchronize();	
	DtoH(d_F2y, F[2]);
	cudaThreadSynchronize();
	DtoH(d_F3y, F[3]);
	cudaThreadSynchronize();
	DtoH(d_F4y, F[4]);
	cudaThreadSynchronize();

	For(i, nt)
	{	
		For(j, nt)
		{	
			For(k, nt)
				cout<<F[1][i*nt*nt + j*nt + k]<<" ";
			
			cout<<endl;
		}
		cout<<endl;
	}
*/

/**************************************************************************************************************************
 * -----------------------------------------------------------------------------------------------------------------------*/      

	derivsZ_Kernel<<<dimGrid, dimBlock>>>(d_W0, d_W0e, d_DW0x, d_DW0y, d_DW0z);
	derivsZ_Kernel<<<dimGrid, dimBlock>>>(d_W1, d_W1e, d_DW1x, d_DW1y, d_DW1z);
	derivsZ_Kernel<<<dimGrid, dimBlock>>>(d_W2, d_W2e, d_DW2x, d_DW2y, d_DW2z);
	derivsZ_Kernel<<<dimGrid, dimBlock>>>(d_W3, d_W3e, d_DW3x, d_DW3y, d_DW3z);
	derivsZ_Kernel<<<dimGrid, dimBlock>>>(d_W4, d_W4e, d_DW4x, d_DW4y, d_DW4z);
	cudaThreadSynchronize();

	flip_Kernel<<<dimGrid, dimBlock>>>(d_W1e, d_W3e, d_DW1x, d_DW3x, d_DW1y, d_DW3y, d_DW1z, d_DW3z);
	cudaThreadSynchronize();
	
	flux<<<dimGrid, dimBlock>>>(d_W0e, d_W1e, d_W2e, d_W3e, d_W4e, d_DW0x, d_DW1x, d_DW2x, d_DW3x, d_DW4x, d_DW0y,
															d_DW1y, d_DW2y, d_DW3y, d_DW4y, d_DW0z, d_DW1z, d_DW2z,
															d_DW3z, d_DW4z, d_F0z, d_F1z, d_F2z, d_F3z, d_F4z, d_tau);
        cudaThreadSynchronize();

	DtoH(d_F0z, F[0]);
	cudaThreadSynchronize();	
	DtoH(d_F1z, F[1]); 
	cudaThreadSynchronize();	
	DtoH(d_F2z, F[2]);
	cudaThreadSynchronize();
	DtoH(d_F3z, F[3]);
	cudaThreadSynchronize();
	DtoH(d_F4z, F[4]);
	cudaThreadSynchronize();

/*	FOR(i, nt)
	{	
		FOR(j, nt)
		{	
			FOR(k, nt)
				cout<<F[2][i*nt*nt + j*nt + k]<<" ";
			
			cout<<endl;
		}
		cout<<endl;
	}
*/

   flipBack_Kernel<<<dimGrid, dimBlock>>>(d_F1z, d_F3z);
   cudaThreadSynchronize();
   
/**********************
 * CHECK with CPU code
 **********************/   
/*	DtoH(d_F0z, F[0]);
	cudaThreadSynchronize();	
	DtoH(d_F1z, F[1]); 
	cudaThreadSynchronize();	
	DtoH(d_F2z, F[2]);
	cudaThreadSynchronize();
	DtoH(d_F3z, F[3]);
	cudaThreadSynchronize();
	DtoH(d_F4z, F[4]);
	cudaThreadSynchronize();

	For(i, nt)
	{	
		For(j, nt)
		{	
			For(k, nt)
				cout<<F[4][i*nt*nt + j*nt + k]<<" ";
			
			cout<<endl;
		}
		cout<<endl;
	}
*/
/**************************************************************************************************************************
 * -----------------------------------------------------------------------------------------------------------------------*/      

/*	update_Kernel<<<dimGrid, dimBlock>>>(d_W0, d_W1, d_W2, d_W3, d_W4, d_F0x, d_F1x, d_F2x, d_F3x, d_F4x,
	                                                                   d_F0y, d_F1y, d_F2y, d_F3y, d_F4y,
	                                                                   d_F0z, d_F1z, d_F2z, d_F3z, d_F4z);
*/
        int block = (nt-1)/(Bsize-1) + 1;

	update_Kernel<<<dim3(block,block,block) , dim3(Bsize, Bsize, Bsize)>>>(d_W0, d_F0x, d_F0y, d_F0z);
	update_Kernel<<<dim3(block,block,block) , dim3(Bsize, Bsize, Bsize)>>>(d_W1, d_F1x, d_F1y, d_F1z);
	update_Kernel<<<dim3(block,block,block) , dim3(Bsize, Bsize, Bsize)>>>(d_W2, d_F2x, d_F2y, d_F2z);
	update_Kernel<<<dim3(block,block,block) , dim3(Bsize, Bsize, Bsize)>>>(d_W3, d_F3x, d_F3y, d_F3z);
	update_Kernel<<<dim3(block,block,block) , dim3(Bsize, Bsize, Bsize)>>>(d_W4, d_F4x, d_F4y, d_F4z);
	
	
	cudaThreadSynchronize();


//	PeriodicBC_Kernel<<<dimGrid, dimBlock>>>(d_W0, d_W1, d_W2, d_W3, d_W4);
	
/** Device to Host Memcpy**/		
/*	DtoH(d_W0, macro->W[0]);	
	DtoH(d_W1, macro->W[1]); 	
	DtoH(d_W2, macro->W[2]);
	DtoH(d_W3, macro->W[3]);
	DtoH(d_W4, macro->W[4]);
	cudaThreadSynchronize();
	
	macro->PeriodicBC();
*/	
	PeriodicBC_Kernel<<<dimGrid, dimBlock>>>(d_W0, d_W1, d_W2, d_W3, d_W4);	

/*	For(i, nt)
	{	
		For(j, nt)
		{	
			For(k, nt)
				cout<<macro->W[4][i*nt*nt + j*nt + k]<<" ";
			cout<<endl;
		}
		cout<<endl;
	}
*/

macro->t = macro->t + macro->dt; 

/***************************************************************************************************************************
 * Kintetic Energy Calculation
 * ************************************************************************************************************************/
	ptype TKE = 0;

	int blks = (Nc-1)/Bsz3D + 1;
	c2p3D_kernel<<<dimGrid, dimBlock>>> (d_W0, d_W1, d_W2, d_W3, d_W4, d_u, d_v, d_w, d_p, d_r, d_k);
	cudaThreadSynchronize();

	KE_kernel<<<blks, Bsz3D>>>(d_k, d_k_reduced);
	cudaThreadSynchronize();
	
	cudaMemcpy(k_reduced, d_k_reduced, sizeof(ptype)*sz_k, cudaMemcpyDeviceToHost);
	
	cout.precision(8);
	
	TKE = 0;	
	FOR(i, sz_k)
	{
		TKE += k_reduced[i];
	}
	
	TKE /= pow(macro->c0*macro->Mt,2);
	TKE /= Nc;

	
/**************************************************************************************************/		
	KE_kernel<<<blks, Bsz3D>>>(d_r, d_k_reduced);
	cudaThreadSynchronize();

	cudaMemcpy(k_reduced, d_k_reduced, sizeof(ptype)*sz_k, cudaMemcpyDeviceToHost);
	cudaThreadSynchronize();
	
	ptype Pmean = 0;
		
	FOR(i, sz_k)
	{
		Pmean += k_reduced[i];
	}
	
	Pmean  /= Nc;
	
	cudaMemcpy(d_Pmean, &Pmean, sizeof(ptype), cudaMemcpyHostToDevice);
	
	PfluctuationSqr_Kernel<<<blks, Bsz3D>>>(d_PflucSqr, d_r, d_Pmean);

	KE_kernel<<<blks, Bsz3D>>>(d_PflucSqr, d_k_reduced);
	cudaThreadSynchronize();

	cudaMemcpy(k_reduced, d_k_reduced, sizeof(ptype)*sz_k, cudaMemcpyDeviceToHost);
	cudaThreadSynchronize();
		
	ptype Prms = 0;
		
	FOR(i, sz_k)
	{
		Prms += k_reduced[i];
	}
	
	Prms = sqrt(Prms/Nc) / (macro->Mt*macro->Mt);

	std::cout<<macro->t/t0<<"   "<<TKE<<"    "<<Prms<< endl;

	
	if(macro->t > 1.56*t0 && m == 0)
	{
		DtoH(d_W0, macro->W[0]);	
		DtoH(d_W1, macro->W[1]); 	
		DtoH(d_W2, macro->W[2]);
		DtoH(d_W3, macro->W[3]);
		DtoH(d_W4, macro->W[4]);
		cudaThreadSynchronize();
		
		macro->SaveVelocityField();
		
		m = m+1;
	}
	
	
	
/**************************************************************************************************************************/


//PrintKE(macro);

}


	DtoH(d_W0, macro->W[0]);	
	DtoH(d_W1, macro->W[1]); 	
	DtoH(d_W2, macro->W[2]);
	DtoH(d_W3, macro->W[3]);
	DtoH(d_W4, macro->W[4]);
	cudaThreadSynchronize();



	CUKILL(d_W0);	CUKILL(d_W1);	CUKILL(d_W2);	CUKILL(d_W3);	CUKILL(d_W4);
	CUKILL(d_W0e);	CUKILL(d_W1e);	CUKILL(d_W2e);	CUKILL(d_W3e);	CUKILL(d_W4e);
		
	CUKILL(d_F0x);	CUKILL(d_F1x);	CUKILL(d_F2x);	CUKILL(d_F3x);	CUKILL(d_F4x);
	CUKILL(d_F0y);	CUKILL(d_F1y);	CUKILL(d_F2y);	CUKILL(d_F3y);	CUKILL(d_F4y);
	CUKILL(d_F0z);	CUKILL(d_F1z);	CUKILL(d_F2z);	CUKILL(d_F3z);	CUKILL(d_F4z);
		
	CUKILL(d_DW0x);	CUKILL(d_DW1x);	CUKILL(d_DW2x);	CUKILL(d_DW3x);	CUKILL(d_DW4x);
	CUKILL(d_DW0y);	CUKILL(d_DW1y);	CUKILL(d_DW2y);	CUKILL(d_DW3y);	CUKILL(d_DW4y);
	CUKILL(d_DW0z);	CUKILL(d_DW1z);	CUKILL(d_DW2z);	CUKILL(d_DW3z);	CUKILL(d_DW4z);
	
	CUKILL(d_tau);
	CUKILL_(d_r);	CUKILL_(d_u);	CUKILL_(d_v);	CUKILL_(d_w);	CUKILL_(d_p);	CUKILL_(d_k);
	
	cudaFree(d_k_reduced);
	cudaFree(d_PflucSqr);
	delete [] k_reduced;

}


/***************************************************************************************************************************
 * ------------------------------------------------------------------------------------------------------------------------*
 * ------------------------------------------------------------------------------------------------------------------------*/

#define Id(i, j, k)  (i)*d_nt*d_nt + (j)*d_nt + (k)
#define IC(i, j, k)  (i)*d_nc*d_nc + (j)*d_nc + (k)

/*
__global__ void derivsX_Kernel(ptype *d_W, ptype *d_We, ptype *d_DWx, ptype *d_DWy, ptype *d_DWz)
{
	int ix = blockIdx.x*blockDim.x + threadIdx.x;
	int iy = blockIdx.y*blockDim.y + threadIdx.y;
	int iz = blockIdx.z*blockDim.z + threadIdx.z;
	 	 	
	ptype WeN, WeS, WeB, WeF;
		
	if(ix < d_nt-1 && iy < d_nt-1 && iz < d_nt-1 && iy > 0 && iz > 0)
	{
		d_We[Id(iy,ix,iz)] = .5*(d_W[Id(iy,ix+1,iz)] + d_W[Id(iy,ix,iz)]);
		
		WeN = .5*(d_W[Id(iy+1,ix+1,iz)] + d_W[Id(iy+1,ix,iz)]);
		
		WeS = .5*(d_W[Id(iy-1,ix+1,iz)] + d_W[Id(iy-1,ix,iz)]);
		
		WeF = .5*(d_W[Id(iy,ix+1,iz+1)] + d_W[Id(iy,ix,iz+1)]);
		
		WeB = .5*(d_W[Id(iy,ix+1,iz-1)] + d_W[Id(iy,ix,iz-1)]);
		
		d_DWx[Id(iy,ix,iz)] = (d_W[Id(iy,ix+1,iz)] - d_W[Id(iy,ix,iz)])/dx;

		__syncthreads();
		
		d_DWy[Id(iy,ix,iz)] = 0.5*(WeN - WeS)/dx;
		
		d_DWz[Id(iy,ix,iz)] = 0.5*(WeF - WeB)/dx;

	}		
	__syncthreads();

}

*/

__global__ void derivsX_Kernel(ptype *d_W, ptype *d_We, ptype *d_DWx, ptype *d_DWy, ptype *d_DWz)
{ 
	int tx = threadIdx.x; 
	int ty = threadIdx.y; 
	int tz = threadIdx.z; 
	
	int ix = blockIdx.x*(Bsize-1) + threadIdx.x;  //try Bsize-2 too !!!
	int iy = blockIdx.y*(Bsize-2) + threadIdx.y; 
	int iz = blockIdx.z*(Bsize-2) + threadIdx.z;
		
	if (ix > d_nt-1 || iy > d_nt-1 || iz > d_nt-1){return;}

	int I = iy*d_nt*d_nt + ix*d_nt + iz;
	
	__shared__ ptype sh_W[Bsize][Bsize][Bsize];
	
	sh_W[ty][tx][tz] = d_W[I];
	__syncthreads();		 
	 	 	
	ptype WeN, WeS, WeB, WeF;
	
	bool bound_check = ( ix < d_nt-1 && iy < d_nt-1 && iz < d_nt-1 && iy > 0 && iz > 0 );
	bool block_check = ( (tx < Bsize-1) && (ty > 0) &&
								(ty < Bsize-1) && (tz > 0) && (tz < Bsize-1) );
			
	if( bound_check && block_check)
	{
		d_We[Id(iy,ix,iz)] = .5*(sh_W[ty][tx+1][tz] + sh_W[ty][tx][tz]);
		
		WeN = .5*(sh_W[ty+1][tx+1][tz] + sh_W[ty+1][tx][tz]);
		
		WeS = .5*(sh_W[ty-1][tx+1][tz] + sh_W[ty-1][tx][tz]);
		
		WeF = .5*(sh_W[ty][tx+1][tz+1] + sh_W[ty][tx][tz+1]);
		
		WeB = .5*(sh_W[ty][tx+1][tz-1] + sh_W[ty][tx][tz-1]);
		
		d_DWx[Id(iy,ix,iz)] = (sh_W[ty][tx+1][tz] - sh_W[ty][tx][tz])/dx;

		__syncthreads();
		
		d_DWy[Id(iy,ix,iz)] = 0.5*(WeN - WeS)/dx;
		
		d_DWz[Id(iy,ix,iz)] = 0.5*(WeF - WeB)/dx;

	}		
	__syncthreads();

}


/*

__global__ void derivsY_Kernel(ptype *d_W, ptype *d_We, ptype *d_DWx, ptype *d_DWy, ptype *d_DWz)
{
	int ix = blockIdx.x*blockDim.x + threadIdx.x;
	int iy = blockIdx.y*blockDim.y + threadIdx.y;
	int iz = blockIdx.z*blockDim.z + threadIdx.z;
	 	 	
	ptype WeN, WeS, WeB, WeF;
		
	if(ix < d_nt-1 && iy < d_nt-1 && iz < d_nt-1 && ix > 0 && iz > 0)
	{
		d_We[Id(iy,ix,iz)] = .5*(d_W[Id(iy+1,ix,iz)] + d_W[Id(iy,ix,iz)]);
		
		WeN = .5*(d_W[Id(iy+1,ix-1,iz)] + d_W[Id(iy,ix-1,iz)]);
		
		WeS = .5*(d_W[Id(iy+1,ix+1,iz)] + d_W[Id(iy,ix+1,iz)]);
		
		WeF = .5*(d_W[Id(iy+1,ix,iz+1)] + d_W[Id(iy,ix,iz+1)]);
		
		WeB = .5*(d_W[Id(iy+1,ix,iz-1)] + d_W[Id(iy,ix,iz-1)]);
		
		d_DWx[Id(iy,ix,iz)] = (d_W[Id(iy+1,ix,iz)] - d_W[Id(iy,ix,iz)])/dx;
		
		__syncthreads();
		
		d_DWy[Id(iy,ix,iz)] = 0.5*(WeS - WeN)/dx;
		
		d_DWz[Id(iy,ix,iz)] = 0.5*(WeF - WeB)/dx;			
	
	}		
	__syncthreads();
	
}

*/


__global__ void derivsY_Kernel(ptype *d_W, ptype *d_We, ptype *d_DWx, ptype *d_DWy, ptype *d_DWz)
{ 
	int tx = threadIdx.x; 
	int ty = threadIdx.y; 
	int tz = threadIdx.z; 
	
	int ix = blockIdx.x*(Bsize-2) + threadIdx.x;  
	int iy = blockIdx.y*(Bsize-2) + threadIdx.y; 
	int iz = blockIdx.z*(Bsize-2) + threadIdx.z;
		
	if (ix > d_nt-1 || iy > d_nt-1 || iz > d_nt-1){return;}

	int I = iy*d_nt*d_nt + ix*d_nt + iz;
	
	__shared__ ptype sh_W[Bsize][Bsize][Bsize];
	
	sh_W[ty][tx][tz] = d_W[I];
	__syncthreads();		 
	 	 	
	ptype WeN, WeS, WeB, WeF;
	
	bool bound_check = ( ix < d_nt-1 && iy < d_nt-1 && iz < d_nt-1 && ix > 0 && iz > 0 );
	bool block_check = ( (ty < Bsize-1) && (tx > 0) &&
								(tx < Bsize-1) && (tz > 0) && (tz < Bsize-1) );
			
	if( bound_check && block_check)
	{
		d_We[Id(iy,ix,iz)] = .5*(sh_W[ty+1][tx][tz] + sh_W[ty][tx][tz]);
		
		WeN = .5*(sh_W[ty+1][tx-1][tz] + sh_W[ty][tx-1][tz]);
		
		WeS = .5*(sh_W[ty+1][tx+1][tz] + sh_W[ty][tx+1][tz]);
		
		WeF = .5*(sh_W[ty+1][tx][tz+1] + sh_W[ty][tx][tz+1]);
		
		WeB = .5*(sh_W[ty+1][tx][tz-1] + sh_W[ty][tx][tz-1]);
		
		d_DWx[Id(iy,ix,iz)] = (sh_W[ty+1][tx][tz] - sh_W[ty][tx][tz])/dx;

		__syncthreads();
		
		d_DWy[Id(iy,ix,iz)] = 0.5*(WeS - WeN)/dx;
		
		d_DWz[Id(iy,ix,iz)] = 0.5*(WeF - WeB)/dx;

	}		
	__syncthreads();

}


/*
__global__ void derivsZ_Kernel(ptype *d_W, ptype *d_We, ptype *d_DWx, ptype *d_DWy, ptype *d_DWz)
{
	int ix = blockIdx.x*blockDim.x + threadIdx.x;
	int iy = blockIdx.y*blockDim.y + threadIdx.y;
	int iz = blockIdx.z*blockDim.z + threadIdx.z;
	 	 	
	ptype WeN, WeS, WeB, WeF;
		
	if(ix < d_nt-1 && iy < d_nt-1 && iz < d_nt-1 && ix > 0 && iy > 0)
	{
		d_We[Id(iy,ix,iz)] = .5*(d_W[Id(iy,ix,iz+1)] + d_W[Id(iy,ix,iz)]);
		
		WeN = .5*(d_W[Id(iy,ix+1,iz+1)] + d_W[Id(iy,ix+1,iz)]);
		
		WeS = .5*(d_W[Id(iy,ix-1,iz+1)] + d_W[Id(iy,ix-1,iz)]);
		
		WeF = .5*(d_W[Id(iy+1,ix,iz+1)] + d_W[Id(iy+1,ix,iz)]);
		
		WeB = .5*(d_W[Id(iy-1,ix,iz+1)] + d_W[Id(iy-1,ix,iz)]);
		
		d_DWx[Id(iy,ix,iz)] = (d_W[Id(iy,ix,iz+1)] - d_W[Id(iy,ix,iz)])/dx;
		
		__syncthreads();
		
		d_DWy[Id(iy,ix,iz)] = (WeN - WeS)/dx;
		
		d_DWz[Id(iy,ix,iz)] = (WeF - WeB)/dx;

	}			
	__syncthreads();

}

*/

__global__ void derivsZ_Kernel(ptype *d_W, ptype *d_We, ptype *d_DWx, ptype *d_DWy, ptype *d_DWz)
{ 
	int tx = threadIdx.x; 
	int ty = threadIdx.y; 
	int tz = threadIdx.z; 
	
	int ix = blockIdx.x*(Bsize-2) + threadIdx.x;  
	int iy = blockIdx.y*(Bsize-2) + threadIdx.y; 
	int iz = blockIdx.z*(Bsize-2) + threadIdx.z; 
		
	if (ix > d_nt-1 || iy > d_nt-1 || iz > d_nt-1){return;}

	int I = iy*d_nt*d_nt + ix*d_nt + iz;
	
	__shared__ ptype sh_W[Bsize][Bsize][Bsize];
	
	sh_W[ty][tx][tz] = d_W[I];
	__syncthreads();		 
	 	 	
	ptype WeN, WeS, WeB, WeF;
	
	bool bound_check = ( ix < d_nt-1 && iy < d_nt-1 && iz < d_nt-1 && ix > 0 && iy > 0 );
	bool block_check = ( (tz < Bsize-1) && (ty > 0) &&
								(ty < Bsize-1) && (tx > 0) && (tx < Bsize-1) );
			
	if( bound_check && block_check)
	{
		d_We[Id(iy,ix,iz)] = .5*(sh_W[ty][tx][tz+1] + sh_W[ty][tx][tz]);
		
		WeN = .5*(sh_W[ty][tx+1][tz+1] + sh_W[ty][tx+1][tz]);
		
		WeS = .5*(sh_W[ty][tx-1][tz+1] + sh_W[ty][tx-1][tz]);
		
		WeF = .5*(sh_W[ty+1][tx][tz+1] + sh_W[ty+1][tx][tz]);
		
		WeB = .5*(sh_W[ty-1][tx][tz+1] + sh_W[ty-1][tx][tz]);
		
		d_DWx[Id(iy,ix,iz)] = (sh_W[ty][tx][tz+1] - sh_W[ty][tx][tz])/dx;

		__syncthreads();
		
		d_DWy[Id(iy,ix,iz)] = 0.5*(WeN - WeS)/dx;
		
		d_DWz[Id(iy,ix,iz)] = 0.5*(WeF - WeB)/dx;

	}		
	__syncthreads();

}


/*
__global__ void flip_Kernel(ptype *d_W1, ptype *d_W2, ptype *d_DW1x, ptype *d_DW2x, ptype *d_DW1y,
                                                      ptype *d_DW2y, ptype *d_DW1z, ptype *d_DW2z)
{
	int ix = blockIdx.x*blockDim.x + threadIdx.x;
	int iy = blockIdx.y*blockDim.y + threadIdx.y;
	int iz = blockIdx.z*blockDim.z + threadIdx.z;
	
	ptype temp;	
	
	int I = iy*d_nt*d_nt + ix*d_nt + iz;
	
	if(ix < d_nt && iy < d_nt && iz < d_nt) {
		
	temp = d_W2[I];			
	__syncthreads();
	
	d_W2[I] = d_W1[I];
	__syncthreads();
	
	d_W1[I] = temp;
	__syncthreads();
	
	temp = d_DW2x[I];			
	__syncthreads();
	
	d_DW2x[I] = d_DW1x[I];	
	__syncthreads();
	
	d_DW1x[I] = temp;			
	__syncthreads();
	
	temp = d_DW2y[I];			
	__syncthreads();
	
	d_DW2y[I] = d_DW1y[I];	
	__syncthreads();
	
	d_DW1y[I] = temp;			
	__syncthreads();
	
	temp = d_DW2z[I];			
	__syncthreads();
	
	d_DW2z[I] = d_DW1z[I];	
	__syncthreads();
	
	d_DW1z[I] = temp;			
	__syncthreads();
	
	}
	
}
*/



__global__ void flip_Kernel(ptype *d_W1, ptype *d_W2, ptype *d_DW1x, ptype *d_DW2x, ptype *d_DW1y,
                                                      ptype *d_DW2y, ptype *d_DW1z, ptype *d_DW2z)
{
	int tx = threadIdx.x; 
	int ty = threadIdx.y; 
	int tz = threadIdx.z; 

	int ix = blockIdx.x*blockDim.x + threadIdx.x;
	int iy = blockIdx.y*blockDim.y + threadIdx.y;
	int iz = blockIdx.z*blockDim.z + threadIdx.z;

	__shared__ ptype sh_W1[Bsize][Bsize][Bsize];
	__shared__ ptype sh_W2[Bsize][Bsize][Bsize];
	__shared__ ptype sh_DW1x[Bsize][Bsize][Bsize];
	__shared__ ptype sh_DW1y[Bsize][Bsize][Bsize];
	__shared__ ptype sh_DW1z[Bsize][Bsize][Bsize];
	__shared__ ptype sh_DW2x[Bsize][Bsize][Bsize];
	__shared__ ptype sh_DW2y[Bsize][Bsize][Bsize];
	__shared__ ptype sh_DW2z[Bsize][Bsize][Bsize];

			
	if(ix < d_nt && iy < d_nt && iz < d_nt) {
	
	int I = iy*d_nt*d_nt + ix*d_nt + iz;
		
	sh_W1[ty][tx][tz] = d_W1[I];			
	sh_W2[ty][tx][tz] = d_W2[I];			
	
	sh_DW1x[ty][tx][tz] = d_DW1x[I];			
	sh_DW2x[ty][tx][tz] = d_DW2x[I];			

	sh_DW1y[ty][tx][tz] = d_DW1y[I];			
	sh_DW2y[ty][tx][tz] = d_DW2y[I];			

	sh_DW1z[ty][tx][tz] = d_DW1z[I];			
	sh_DW2z[ty][tx][tz] = d_DW2z[I];			

	__syncthreads();
	
	d_W1[I] = sh_W2[ty][tx][tz];			
	d_W2[I] = sh_W1[ty][tx][tz];			
	
	d_DW1x[I] = sh_DW2x[ty][tx][tz];			
	d_DW2x[I] = sh_DW1x[ty][tx][tz];			

	d_DW1y[I] = sh_DW2y[ty][tx][tz];			
	d_DW2y[I] = sh_DW1y[ty][tx][tz];			
	
	d_DW1z[I] = sh_DW2z[ty][tx][tz];			
	d_DW2z[I] = sh_DW1z[ty][tx][tz];			
	
	__syncthreads();	
	}
	
}



/*
__global__ void flipBack_Kernel(ptype *d_F1, ptype *d_F2)
{
	int ix = blockIdx.x*blockDim.x + threadIdx.x;
	int iy = blockIdx.y*blockDim.y + threadIdx.y;
	int iz = blockIdx.z*blockDim.z + threadIdx.z;
	
	int I = iy*d_nt*d_nt + ix*d_nt + iz;
	
	if(ix < d_nt && iy < d_nt && iz < d_nt) {
	
	ptype temp;
	
	temp = d_F1[I];			
	__syncthreads();
	
	d_F1[I] = d_F2[I];		
	__syncthreads();
	
	d_F2[I] = temp;			
	__syncthreads();
	
	}
}

*/


__global__ void flipBack_Kernel(ptype *d_F1, ptype *d_F2)
{
	int tx = threadIdx.x; 
	int ty = threadIdx.y; 
	int tz = threadIdx.z; 

	int ix = blockIdx.x*blockDim.x + threadIdx.x;
	int iy = blockIdx.y*blockDim.y + threadIdx.y;
	int iz = blockIdx.z*blockDim.z + threadIdx.z;

	__shared__ ptype sh_F1[Bsize][Bsize][Bsize];
	__shared__ ptype sh_F2[Bsize][Bsize][Bsize];
	

	
	if(ix < d_nt && iy < d_nt && iz < d_nt) {
	
	int I = iy*d_nt*d_nt + ix*d_nt + iz;	
	
	sh_F1[ty][tx][tz] = d_F1[I];		
	sh_F2[ty][tx][tz] = d_F2[I];		

	__syncthreads();
	
	d_F1[I] = sh_F2[ty][tx][tz];			
	d_F2[I] = sh_F1[ty][tx][tz];			

	__syncthreads();
	
	}
}



/*
__global__ void update_Kernel (ptype *d_W0, ptype *d_W1, ptype *d_W2, ptype *d_W3, ptype *d_W4,
						       ptype *d_F0x, ptype *d_F1x, ptype *d_F2x, ptype *d_F3x, ptype *d_F4x,
						       ptype *d_F0y, ptype *d_F1y, ptype *d_F2y, ptype *d_F3y, ptype *d_F4y,
						       ptype *d_F0z, ptype *d_F1z, ptype *d_F2z, ptype *d_F3z, ptype *d_F4z)

{
	int ix = blockIdx.x*blockDim.x + threadIdx.x;
	int iy = blockIdx.y*blockDim.y + threadIdx.y;
	int iz = blockIdx.z*blockDim.z + threadIdx.z;
	int I  = iy*d_nt*d_nt + ix*d_nt + iz;

	if(ix < d_nt-1 && iy < d_nt-1 && iz < d_nt-1 && ix > 0 && iy > 0 && iz > 0)
	{
		
		d_W0[I] = d_W0[I] - (1/dx) * ( d_F0x[Id(iy, ix, iz)] - d_F0x[Id(iy, ix-1, iz)] +
		                               d_F0y[Id(iy, ix, iz)] - d_F0y[Id(iy-1, ix, iz)] +
											    d_F0z[Id(iy, ix, iz)] - d_F0z[Id(iy, ix, iz-1)] );
		

		d_W1[I] = d_W1[I] - (1/dx) * ( d_F1x[Id(iy, ix, iz)] - d_F1x[Id(iy, ix-1, iz)] +
											    d_F1y[Id(iy, ix, iz)] - d_F1y[Id(iy-1, ix, iz)] +
											    d_F1z[Id(iy, ix, iz)] - d_F1z[Id(iy, ix, iz-1)] );
		
		
		d_W2[I] = d_W2[I] - (1/dx) * ( d_F2x[Id(iy, ix, iz)] - d_F2x[Id(iy, ix-1, iz)] +
											    d_F2y[Id(iy, ix, iz)] - d_F2y[Id(iy-1, ix, iz)] +
											    d_F2z[Id(iy, ix, iz)] - d_F2z[Id(iy, ix, iz-1)] );
		

		d_W3[I] = d_W3[I] - (1/dx) * ( d_F3x[Id(iy, ix, iz)] - d_F3x[Id(iy, ix-1, iz)] +
											    d_F3y[Id(iy, ix, iz)] - d_F3y[Id(iy-1, ix, iz)] +
											    d_F3z[Id(iy, ix, iz)] - d_F3z[Id(iy, ix, iz-1)] );

		
		d_W4[I] = d_W4[I] - (1/dx) * ( d_F4x[Id(iy, ix, iz)] - d_F4x[Id(iy, ix-1, iz)] +
												 d_F4y[Id(iy, ix, iz)] - d_F4y[Id(iy-1, ix, iz)] +
												 d_F4z[Id(iy, ix, iz)] - d_F4z[Id(iy, ix, iz-1)] );											
	
	}
	
	__syncthreads();

}

*/

/*
__global__ void update_Kernel (ptype *d_W, ptype *d_Fx, ptype *d_Fy, ptype *d_Fz)
{
	int ix = blockIdx.x*blockDim.x + threadIdx.x;
	int iy = blockIdx.y*blockDim.y + threadIdx.y;
	int iz = blockIdx.z*blockDim.z + threadIdx.z;
	int I  = iy*d_nt*d_nt + ix*d_nt + iz;

	if(ix < d_nt-1 && iy < d_nt-1 && iz < d_nt-1 && ix > 0 && iy > 0 && iz > 0)
	{
		
		d_W[I] = d_W[I] - (1/dx) * ( d_Fx[Id(iy, ix, iz)] - d_Fx[Id(iy, ix-1, iz)] +
		                               d_Fy[Id(iy, ix, iz)] - d_Fy[Id(iy-1, ix, iz)] +
											    d_Fz[Id(iy, ix, iz)] - d_Fz[Id(iy, ix, iz-1)] );
		
	}
	
	__syncthreads();

}

*/




/*
__global__ void update_Kernel (ptype *d_W0, ptype *d_W1, ptype *d_W2, ptype *d_W3, ptype *d_W4,
						       ptype *d_F0x, ptype *d_F1x, ptype *d_F2x, ptype *d_F3x, ptype *d_F4x,
						       ptype *d_F0y, ptype *d_F1y, ptype *d_F2y, ptype *d_F3y, ptype *d_F4y,
						       ptype *d_F0z, ptype *d_F1z, ptype *d_F2z, ptype *d_F3z, ptype *d_F4z)

{
	int tx = threadIdx.x; 
	int ty = threadIdx.y; 
	int tz = threadIdx.z; 
	
	int ix = blockIdx.x*(Bsize-1) + threadIdx.x;
	int iy = blockIdx.y*(Bsize-1) + threadIdx.y; 
	int iz = blockIdx.z*(Bsize-1) + threadIdx.z;
		
	if (ix > d_nt-1 || iy > d_nt-1 || iz > d_nt-1){return;}

	int I = iy*d_nt*d_nt + ix*d_nt + iz;
	
	__shared__ ptype sh_F0x[Bsize][Bsize][Bsize];
	__shared__ ptype sh_F0y[Bsize][Bsize][Bsize];
	__shared__ ptype sh_F0z[Bsize][Bsize][Bsize];
	
	__shared__ ptype sh_F1x[Bsize][Bsize][Bsize];
	__shared__ ptype sh_F1y[Bsize][Bsize][Bsize];
	__shared__ ptype sh_F1z[Bsize][Bsize][Bsize];

	__shared__ ptype sh_F2x[Bsize][Bsize][Bsize];
	__shared__ ptype sh_F2y[Bsize][Bsize][Bsize];
	__shared__ ptype sh_F2z[Bsize][Bsize][Bsize];

	__shared__ ptype sh_F3x[Bsize][Bsize][Bsize];
	__shared__ ptype sh_F3y[Bsize][Bsize][Bsize];
	__shared__ ptype sh_F3z[Bsize][Bsize][Bsize];

	__shared__ ptype sh_F4x[Bsize][Bsize][Bsize];
	__shared__ ptype sh_F4y[Bsize][Bsize][Bsize];
	__shared__ ptype sh_F4z[Bsize][Bsize][Bsize];


	sh_F0x[ty][tx][tz] = d_F0x[I];
	sh_F0y[ty][tx][tz] = d_F0y[I];
	sh_F0z[ty][tx][tz] = d_F0z[I];
	
	sh_F1x[ty][tx][tz] = d_F1x[I];
	sh_F1y[ty][tx][tz] = d_F1y[I];
	sh_F1z[ty][tx][tz] = d_F1z[I];

	sh_F2x[ty][tx][tz] = d_F2x[I];
	sh_F2y[ty][tx][tz] = d_F2y[I];
	sh_F2z[ty][tx][tz] = d_F2z[I];

	sh_F3x[ty][tx][tz] = d_F3x[I];
	sh_F3y[ty][tx][tz] = d_F3y[I];
	sh_F3z[ty][tx][tz] = d_F3z[I];

	sh_F4x[ty][tx][tz] = d_F4x[I];
	sh_F4y[ty][tx][tz] = d_F4y[I];
	sh_F4z[ty][tx][tz] = d_F4z[I];
	
	
	__syncthreads();		 
	 	 	
	 	 		
	bool bound_check = ( ix < d_nt-1 && iy < d_nt-1 && iz < d_nt-1 && ix > 0 && iy > 0 && iz > 0 );
	bool block_check = ( (tx > 0) && (tx < Bsize-1) && (ty > 0) );
			
	if( bound_check && block_check)	
	{
		I = iy*d_nt*d_nt + ix*d_nt + iz;
		
		d_W0[I] = d_W0[I] - (1/dx) * ( sh_F0x[ty][tx][tz] - sh_F0x[ty][tx-1][tz] +
		                               sh_F0y[ty][tx][tz] - sh_F0y[ty-1][tx][tz] +
											    sh_F0z[ty][tx][tz] - sh_F0z[ty][tx][tz-1] );
	
		
		d_W1[I] = d_W1[I] - (1/dx) * ( sh_F1x[ty][tx][tz] - sh_F1x[ty][tx-1][tz] +
		                               sh_F1y[ty][tx][tz] - sh_F1y[ty-1][tx][tz] +
											    sh_F1z[ty][tx][tz] - sh_F1z[ty][tx][tz-1] );


		d_W2[I] = d_W2[I] - (1/dx) * ( sh_F2x[ty][tx][tz] - sh_F2x[ty][tx-1][tz] +
		                               sh_F2y[ty][tx][tz] - sh_F2y[ty-1][tx][tz] +
											    sh_F2z[ty][tx][tz] - sh_F2z[ty][tx][tz-1] );


		d_W3[I] = d_W3[I] - (1/dx) * ( sh_F3x[ty][tx][tz] - sh_F3x[ty][tx-1][tz] +
		                               sh_F3y[ty][tx][tz] - sh_F3y[ty-1][tx][tz] +
											    sh_F3z[ty][tx][tz] - sh_F3z[ty][tx][tz-1] );


		d_W4[I] = d_W4[I] - (1/dx) * ( sh_F4x[ty][tx][tz] - sh_F4x[ty][tx-1][tz] +
		                               sh_F4y[ty][tx][tz] - sh_F4y[ty-1][tx][tz] +
											    sh_F4z[ty][tx][tz] - sh_F4z[ty][tx][tz-1] );

	}
	
	__syncthreads();

}

*/


__global__ void update_Kernel (ptype *d_W, ptype *d_Fx, ptype *d_Fy, ptype *d_Fz)
{
	int tx = threadIdx.x; 
	int ty = threadIdx.y; 
	int tz = threadIdx.z; 
	
	int ix = blockIdx.x*(Bsize-1) + threadIdx.x;
	int iy = blockIdx.y*(Bsize-1) + threadIdx.y; 
	int iz = blockIdx.z*(Bsize-1) + threadIdx.z;
		
	if (ix > d_nt-1 || iy > d_nt-1 || iz > d_nt-1)	{return;}

	int I = iy*d_nt*d_nt + ix*d_nt + iz;
	
	__shared__ ptype sh_Fx[Bsize][Bsize][Bsize];
	__shared__ ptype sh_Fy[Bsize][Bsize][Bsize];
	__shared__ ptype sh_Fz[Bsize][Bsize][Bsize];
	__shared__ ptype sh_W [Bsize][Bsize][Bsize];
	

	sh_Fx[ty][tx][tz] = d_Fx[I];
	sh_Fy[ty][tx][tz] = d_Fy[I];
	sh_Fz[ty][tx][tz] = d_Fz[I];
	
	sh_W[ty][tx][tz]  = d_W[I];


	__syncthreads();		 
	 	 	
	 	 		
	bool bound_check = ( ix < d_nt-1 && iy < d_nt-1 && iz < d_nt-1 && ix > 0 && iy > 0 && iz > 0 );
	bool block_check = ( (tx > 0) && (tz > 0) && (ty > 0) );
			
	if( bound_check && block_check)	
	{
		I = iy*d_nt*d_nt + ix*d_nt + iz;
		
		sh_W[ty][tx][tz] = sh_W[ty][tx][tz] - (1/dx) * ( sh_Fx[ty][tx][tz] - sh_Fx[ty][tx-1][tz] +
		                             sh_Fy[ty][tx][tz] - sh_Fy[ty-1][tx][tz] +
											  sh_Fz[ty][tx][tz] - sh_Fz[ty][tx][tz-1] );
											  	
		__syncthreads();
											  
		d_W[I] = sh_W[ty][tx][tz];
				
	}
	
	__syncthreads();
	
}


__global__ void PeriodicBC_Kernel (ptype *d_W0, ptype *d_W1, ptype *d_W2, ptype *d_W3, ptype *d_W4)
{
	int ix = blockIdx.x*blockDim.x + threadIdx.x;
	int iy = blockIdx.y*blockDim.y + threadIdx.y;
	int iz = blockIdx.z*blockDim.z + threadIdx.z;
	
	
   if(ix < d_nt && iy < d_nt && iz < d_nt) {
	
	d_W0[Id(iy, 0, iz)]       =  d_W0[Id(iy, d_nt-2, iz)];
	d_W0[Id(iy, d_nt-1, iz)]  =  d_W0[Id(iy, 1, iz)];
	
	d_W1[Id(iy, 0, iz)]       =  d_W1[Id(iy, d_nt-2, iz)];
	d_W1[Id(iy, d_nt-1, iz)]  =  d_W1[Id(iy, 1, iz)];

	d_W2[Id(iy, 0, iz)]       =  d_W2[Id(iy, d_nt-2, iz)];
	d_W2[Id(iy, d_nt-1, iz)]  =  d_W2[Id(iy, 1, iz)];

	d_W3[Id(iy, 0, iz)] 		  =  d_W3[Id(iy, d_nt-2, iz)];
	d_W3[Id(iy, d_nt-1, iz)]  =  d_W3[Id(iy, 1, iz)];

	d_W4[Id(iy, 0, iz)]       =  d_W4[Id(iy, d_nt-2, iz)];
	d_W4[Id(iy, d_nt-1, iz)]  =  d_W4[Id(iy, 1, iz)];


//
	d_W0[Id(0, ix, iz)]       =  d_W0[Id(d_nt-2, ix, iz)];
	d_W0[Id(d_nt-1, ix, iz)]  =  d_W0[Id(1, ix, iz)];
	
	d_W1[Id(0, ix, iz)]       =  d_W1[Id(d_nt-2, ix, iz)];
	d_W1[Id(d_nt-1, ix, iz)]  =  d_W1[Id(1, ix, iz)];

	d_W2[Id(0, ix, iz)]       =  d_W2[Id(d_nt-2, ix, iz)];
	d_W2[Id(d_nt-1, ix, iz)]  =  d_W2[Id(1, ix, iz)];

	d_W3[Id(0, ix, iz)]       =  d_W3[Id(d_nt-2, ix, iz)];
	d_W3[Id(d_nt-1, ix, iz)]  =  d_W3[Id(1, ix, iz)];

	d_W4[Id(0, ix, iz)]       =  d_W4[Id(d_nt-2, ix, iz)];
	d_W4[Id(d_nt-1, ix, iz)]  =  d_W4[Id(1, ix, iz)];
	
//
	d_W0[Id(iy, ix, 0)]       =  d_W0[Id(iy, ix, d_nt-2)];
	d_W0[Id(iy, ix, d_nt-1)]  =  d_W0[Id(iy, ix, 1)];

	d_W1[Id(iy, ix, 0)]       =  d_W1[Id(iy, ix, d_nt-2)];
	d_W1[Id(iy, ix, d_nt-1)]  =  d_W1[Id(iy, ix, 1)];

	d_W2[Id(iy, ix, 0)]       =  d_W2[Id(iy, ix, d_nt-2)];
	d_W2[Id(iy, ix, d_nt-1)]  =  d_W2[Id(iy, ix, 1)];

	d_W3[Id(iy, ix, 0)]       =  d_W3[Id(iy, ix, d_nt-2)];
	d_W3[Id(iy, ix, d_nt-1)]  =  d_W3[Id(iy, ix, 1)];

	d_W4[Id(iy, ix, 0)]       =  d_W4[Id(iy, ix, d_nt-2)];
	d_W4[Id(iy, ix, d_nt-1)]  =  d_W4[Id(iy, ix, 1)];

	}
	__syncthreads();
}



/***************************************************************************************************************************
 * ------------------------------------------------------------------------------------------------------------------------*/
__global__ void flux(ptype *d_W0, ptype *d_W1, ptype *d_W2, ptype *d_W3, ptype *d_W4,
                     ptype *d_DW0x, ptype *d_DW1x, ptype *d_DW2x, ptype *d_DW3x, ptype *d_DW4x,
							ptype *d_DW0y, ptype *d_DW1y, ptype *d_DW2y, ptype *d_DW3y, ptype *d_DW4y,
							ptype *d_DW0z, ptype *d_DW1z, ptype *d_DW2z, ptype *d_DW3z, ptype *d_DW4z,
							ptype *d_F0, ptype *d_F1, ptype *d_F2, ptype *d_F3, ptype *d_F4, ptype *d_tau)

{
	
	int ix = blockIdx.x*blockDim.x + threadIdx.x;
	int iy = blockIdx.y*blockDim.y + threadIdx.y;
	int iz = blockIdx.z*blockDim.z + threadIdx.z;
	
	if(ix < d_nt && iy < d_nt && iz < d_nt) {
					
	int Id = (iy)*d_nt*d_nt + (ix)*d_nt + (iz);


	ptype r, u, v, w, p, T, lam, U[3];
	
	r = d_W0[Id];
	__syncthreads();
	
	u = d_W1[Id] / r;

	v = d_W2[Id] / r;

	w = d_W3[Id] / r;
	
	__syncthreads();
	
	T = (GAM-1) * ( d_W4[Id]/r - .5*(u*u + v*v +w*w) );	
	__syncthreads();
	

	p = T * r;
	__syncthreads();
	

	lam = .5*r / p;

	
	U[0] = u; U[1] = v; U[2] = w;
	__syncthreads();
		

//	ptype tau = mu0 * (pow((T/T0),1.5)*(T0+110.4)/(T+110.4)) / p; 
	ptype tau = mu0 * (pow((double(T/T0)),0.76)) / p; 

/* -------------------------------------------------------------------*
 * Moment Integrals Calculation
 * -------------------------------------------------------------------*/
	
	ptype Ie20,Ie40; //Ie00
	ptype If0[3*7];

//	Ie00 = 1;
	Ie20 = K/(2*lam);
	Ie40 = 3*K/(4*lam*lam) + K*(K-1)/(4*lam*lam);
	
	__syncthreads();
	
	for (int j=0;j<3;j++)
	{
		If0[j*7 + 0] = 1;
		If0[j*7 + 1] = U[j];
		__syncthreads();
		
		for (int i=2;i<7;i++)
		{
			If0[j*7 + i] = U[j]*If0[j*7 + i-1] + If0[j*7 + i-2]*(i-1)/(2*lam);
			
			__syncthreads();
		}
	}		
	
/* ---------------------------------------------------------*
 * MOMENT MATRICES CALCULATION
 * ---------------------------------------------------------*/

	ptype Mu_ax[5];	ptype Mv_ay[5];	ptype Mw_az[5];
	ptype Mu_A [5];	ptype M[5];
		
	ptype Muv_ay[5];	ptype Muw_az[5];	ptype Mu2_ax[5];	
	
	__syncthreads();
		
/* -------------------------------------------------*
 * Slope calculation
 * -------------------------------------------------*/	
	
	ptype bx[5], by[5], bz[5], ax[5], ay[5], az[5];
	ptype A[5],  B[5];
				
	bx[0] = d_DW0x[Id] / r;
	by[0] = d_DW0y[Id] / r;
	bz[0] = d_DW0z[Id] / r;

	bx[1] = d_DW1x[Id] / r;
	by[1] = d_DW1y[Id] / r;
	bz[1] = d_DW1z[Id] / r;

	bx[2] = d_DW2x[Id] / r;
	by[2] = d_DW2y[Id] / r;
	bz[2] = d_DW2z[Id] / r;

	bx[3] = d_DW3x[Id] / r;
	by[3] = d_DW3y[Id] / r;
	bz[3] = d_DW3z[Id] / r;

	bx[4] = d_DW4x[Id] / r;
	by[4] = d_DW4y[Id] / r;
	bz[4] = d_DW4z[Id] / r;
	
	__syncthreads();
	
	d_slopesolver(bx, U, lam, ax);
	d_slopesolver(by, U, lam, ay);
	d_slopesolver(bz, U, lam, az);
	
	__syncthreads();

// A Calculation


	d_MCal(Mu_ax,If0,Ie20,Ie40,1,0,0, ax);	
	d_MCal(Mv_ay,If0,Ie20,Ie40,0,1,0, ay);	
	d_MCal(Mw_az,If0,Ie20,Ie40,0,0,1, az);
	
	d_MCal(Muv_ay, If0, Ie20, Ie40,1,1,0, ay);	
	d_MCal(Muw_az, If0, Ie20, Ie40,1,0,1, az);	
	d_MCal(Mu2_ax, If0, Ie20, Ie40,2,0,0, ax);
	
	__syncthreads();
	

	for (int i = 0; i < 5; i++)
	{
		B[i] = -Mu_ax[i] - Mv_ay[i] - Mw_az[i];
	}
	__syncthreads();
	

	d_slopesolver(B, U, lam, A);
	
	__syncthreads();

	d_MCal(Mu_A,If0,Ie20,Ie40,1,0,0, A);	
	
	__syncthreads();
	
	
	ptype val0	= .5 * ( If0[0*7 + 3] + If0[0*7 + 1]*If0[1*7 + 2] + If0[0*7 + 1]*If0[2*7 + 2] + If0[0*7 + 1]*Ie20  );
	
	__syncthreads();
	
	M[0] = If0[0*7 + 1];
	M[1] = If0[0*7 + 2];
	M[2] = If0[0*7 + 1]*If0[1*7 + 1];
	M[3] = If0[0*7 + 1]*If0[2*7 + 1];
	M[4] = val0;	
	
	__syncthreads();

/* ----------------------------------------------*
 * Flux calculation
 * ----------------------------------------------*/

// integral dt
	
	ptype p0, p1, p2;
	
	p0 = dt;	p1 = -tau*dt;	p2 = dt*dt/2;
	__syncthreads();

	
	__syncthreads();
	
	ptype F[5] = {0};
	__syncthreads();
		
	FOR(i, 5) {
		F[i] = r * ( (p1+p2)*Mu_A[i] + p1*Mu2_ax[i] + p1*Muv_ay[i] + p1*Muw_az[i] +  p0*M[i]);  
	}

	__syncthreads();
	
	d_F0[Id] = F[0];

	d_F1[Id] = F[1];

	d_F2[Id] = F[2];

	d_F3[Id] = F[3];

	d_F4[Id] = F[4]; 		


	__syncthreads();
		

	} // end of if statement for checking thread out of range acess

}
	
/***************************************************************************************************************************/	
	
/* ----------------------------------------------------------------------*
 *  SLOPE SOLVER
 * ----------------------------------------------------------------------*/

__device__ void d_slopesolver(ptype b[5], ptype U[3], ptype lam, ptype a[5])
{
	ptype R2, R3, R4, R5;
	
	R2 = b[1] - U[0]*b[0];
	R3 = b[2] - U[1]*b[0];
	R4 = b[3] - U[2]*b[0];
	R5 = 2*b[4] - b[0]*(U[0]*U[0]+U[1]*U[1]+U[2]*U[2]+(K+3)/(2*lam));
	
	__syncthreads();
	
	a[4] = (1/PRN)*(R5-2*U[0]*R2-2*U[1]*R3-2*U[2]*R4)*(4*lam*lam)/(K+3);
	__syncthreads();
	a[3] = 2*lam*R4 - U[2]*a[4];
	__syncthreads();
	a[2] = 2*lam*R3 - U[1]*a[4];
	__syncthreads();
	a[1] = 2*lam*R2 - U[0]*a[4];
	__syncthreads();
	
	__syncthreads();
	
	a[0] = b[0] - a[1]*U[0] - a[2]*U[1] -a[3]*U[2]-.5*a[4]*(U[0]*U[0] +
	       U[1]*U[1] + U[2]*U[2]+(K+3)/(2*lam));
	
	__syncthreads();
}

/* ----------------------------------------------------------------------*
 *  Moment - Matrix - Multiplication
 * ----------------------------------------------------------------------*/

__device__ void d_mult(ptype M[5*2], int k, ptype C, ptype B[5])
{
	for (int i=0;i<5;i++)
	{
		B[i]=C*M[i*2+k];
	}
	__syncthreads();
}

__device__ void d_Mult(ptype M[5*5], ptype C, ptype A[5], ptype B[5])
{
	B[0]=0;B[1]=0;B[2]=0;B[3]=0;B[4]=0;
	
	for (int i=0;i<5;i++)
	{
		__syncthreads();
		
		B[i] = B[i] + C * M[i*5 + i] * A[i];
		
		for (int j=i+1;j<5;j++)
		{
			__syncthreads();
			
			B[i] = B[i] + C * M[i*5 + j] * A[j];
			__syncthreads();
			
			B[j] = B[j] + C * M[i*5 + j] * A[i];
		}
	}
	__syncthreads();
}


// Vector Addition function
__device__ void d_add(ptype A[5], ptype B[5])
{
	for (int i=0;i<5;i++)
	{
		B[i] = A[i] + B[i];
		__syncthreads();
	}
}

/* -------------------------------------------------------------------------------------*
 * 5x5 Moment calculator
 * -------------------------------------------------------------------------------------*/

__device__ void d_MCal(ptype M[5],ptype If[3*7],ptype Ie2,ptype Ie4, int k,int l,int m, ptype ax[5])
{
	ptype val0, val1, val2, val3;
	 
	val0 = .5 * ( If[0*7 + 2+k]*If[1*7 + l]*If[2*7 + m] + If[0*7 + k]*If[1*7 + 2+l]*If[2*7 + m] +
	              If[0*7 + k]*If[1*7 + l]*If[2*7 + 2+m] + If[0*7 + k]*If[1*7 + l]*If[2*7 + m]*Ie2  );

	__syncthreads();
	
	M[0] = ax[0]*(If[0*7 + k]*If[1*7 + l]*If[2*7 + m]) + ax[1]*(If[0*7 + 1+k]*If[1*7 + l]*If[2*7 + m]) +
			 ax[2]*(If[0*7 + k]*If[1*7 + 1+l]*If[2*7 + m]) +
	       ax[3]*(If[0*7 + k]*If[1*7 + l]*If[2*7 + 1+m]) + ax[4]*val0 ;	


				  
	val1 = .5 * ( If[0*7 + 3+k]*If[1*7 + l]*If[2*7 + m] + If[0*7 + 1+k]*If[1*7 + 2+l]*If[2*7 + m] +
	              If[0*7 + 1+k]*If[1*7 + l]*If[2*7 + 2+m] +
					  If[0*7 + 1+k]*If[1*7 + l]*If[2*7 + m]*Ie2  );
	

	__syncthreads();

	
	M[1] = ax[0]*(If[0*7 + 1+k]*If[1*7 + l]*If[2*7 + m]) + ax[1]*(If[0*7 + 2+k]*If[1*7 + l]*If[2*7 + m]) +
	       ax[2]*(If[0*7 + 1+k]*If[1*7 + 1+l]*If[2*7 + m]) +
			 ax[3]*(If[0*7 + 1+k]*If[1*7 + l]*If[2*7 + 1+m]) + ax[4]*val1;
				 

   
   val2 = .5 * ( If[0*7 + 2+k]*If[1*7 + 1+l]*If[2*7 + m] + If[0*7 + k]*If[1*7 + 3+l]*If[2*7 + m] + 
                 If[0*7 + k]*If[1*7 + 1+l]*If[2*7 + 2+m] +
					  If[0*7 + k]*If[1*7 + 1+l]*If[2*7 + m]*Ie2  );
	
	__syncthreads();
	
	
	M[2] = ax[0]*(If[0*7 + k]*If[1*7 + 1+l]*If[2*7 + m])   + ax[1]*(If[0*7 + 1+k]*If[1*7 + 1+l]*If[2*7 + m]) +
	       ax[2]*(If[0*7 + k]*If[1*7 + 2+l]*If[2*7 + m]) +
			 ax[3]*(If[0*7 + k]*If[1*7 + 1+l]*If[2*7 + 1+m]) + ax[4]*val2;
				 
	

	val3 = .5 * ( If[0*7 + 2+k]*If[1*7 + l]*If[2*7 + 1+m] + If[0*7 + k]*If[1*7 + 2+l]*If[2*7 + 1+m] +
	              If[0*7 + k]*If[1*7 + l]*If[2*7 + 3+m] + 
	              If[0*7 + k]*If[1*7 + l]*If[2*7 + 1+m]*Ie2  );
	
	__syncthreads();
	
	M[3] = ax[0]*(If[0*7 + k]*If[1*7 + l]*If[2*7 + 1+m]) + ax[1]*(If[0*7 + 1+k]*If[1*7 + l]*If[2*7 + 1+m]) +
	       ax[2]*(If[0*7 + k]*If[1*7 + 1+l]*If[2*7 + 1+m]) +
	       ax[3]*(If[0*7 + k]*If[1*7 + l]*If[2*7 + 2+m]) + ax[4]*val3;
				 


	M[4] = .25*ax[4]* ( If[0*7 + 4+k]*If[1*7 + l]*If[2*7 + m] + If[0*7 + k]*If[1*7 + 4+l]*If[2*7 + m] +
	                    If[0*7 + k]*If[1*7 + l]*If[2*7 + 4+m] + If[0*7 + k]*If[1*7 + l]*If[2*7 + m]*Ie4 +
	                  2*If[0*7 + 2+k]*If[1*7 + 2+l]*If[2*7 + m] + 2*If[0*7 + 2+k]*If[1*7 + l]*If[2*7 + 2+m] +
	                  2*If[0*7 + 2+k]*If[1*7 + l]*If[2*7 + m]*Ie2 + 2*If[0*7 + k]*If[1*7 + 2+l]*If[2*7 + 2+m] +
	                  2*If[0*7 + k]*If[1*7 + 2+l]*If[2*7 + m]*Ie2 + 2*If[0*7 + k]*If[1*7 + l]*If[2*7 + 2+m]*Ie2 ) +
			 ax[0]*val0 + ax[1]*val1 + ax[2]*val2 + ax[3]*val3;


	__syncthreads();
}



void PrintKE(Macro *macro)
{
	std::cout.precision(10);
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


__global__ void c2p3D_kernel (ptype *W0,ptype *W1,ptype *W2,ptype *W3,ptype *W4,
										ptype *u, ptype *v, ptype *w, ptype *p, ptype *r, ptype *k)
{
	int tx = threadIdx.x;
	int ty = threadIdx.y;
	int tz = threadIdx.z;
	
	int ix = blockIdx.x*blockDim.x + threadIdx.x;
	int iy = blockIdx.y*blockDim.y + threadIdx.y;
	int iz = blockIdx.z*blockDim.z + threadIdx.z;

	__shared__ ptype R[Bsize][Bsize][Bsize];
	__shared__ ptype U[Bsize][Bsize][Bsize];
	__shared__ ptype V[Bsize][Bsize][Bsize];
	__shared__ ptype W[Bsize][Bsize][Bsize];
	__shared__ ptype E[Bsize][Bsize][Bsize];
	__shared__ ptype KE[Bsize][Bsize][Bsize];
	
	
	if( ix < d_nt-2 && iy < d_nt-2 && iz < d_nt-2) {

	R[ty][tx][tz] = W0[Id(iy+1, ix+1, iz+1)];
	__syncthreads();
	
	U[ty][tx][tz] = W1[Id(iy+1, ix+1, iz+1)] / R[ty][tx][tz];
	V[ty][tx][tz] = W2[Id(iy+1, ix+1, iz+1)] / R[ty][tx][tz];
	W[ty][tx][tz] = W3[Id(iy+1, ix+1, iz+1)] / R[ty][tx][tz];
	E[ty][tx][tz] = W4[Id(iy+1, ix+1, iz+1)];
	
	__syncthreads();
	
	KE[ty][tx][tz] = R[ty][tx][tz]*(U[ty][tx][tz]*U[ty][tx][tz] + V[ty][tx][tz]*V[ty][tx][tz] + W[ty][tx][tz]*W[ty][tx][tz]);

	}	
	

	__syncthreads();
			
	if(ix < d_nt-2 && iy < d_nt-2 && iz < d_nt-2)
	{
		r[IC(iy, ix, iz)]  =  R[ty][tx][tz];
		
		u[IC(iy, ix, iz)]  =  U[ty][tx][tz];
		v[IC(iy, ix, iz)]  =  V[ty][tx][tz];
		w[IC(iy, ix, iz)]  =  W[ty][tx][tz];
		k[IC(iy, ix, iz)]  =  KE[ty][tx][tz];
		
		p[IC(iy, ix, iz)]  =  (GAM-1) * ( E[ty][tx][tz] - 0.5 * R[ty][tx][tz] * KE[ty][tx][tz] ); 
	}
	
	__syncthreads();	
}

/*
__global__ void c2p3D_kernel (ptype *W0,ptype *W1,ptype *W2,ptype *W3,ptype *W4,
										ptype *u, ptype *v, ptype *w, ptype *p, ptype *r, ptype *k)
{
	int ix = blockIdx.x*blockDim.x + threadIdx.x;
	int iy = blockIdx.y*blockDim.y + threadIdx.y;
	int iz = blockIdx.z*blockDim.z + threadIdx.z;

		
	if( ix < d_nt-1 && iy < d_nt-1 && iz < d_nt-1 && ix > 0 && iy > 0 && iz > 0) {

	r[IC(iy-1, ix-1, iz-1)] = W0[Id(iy, ix, iz)];
	__syncthreads();
	
	u[IC(iy-1, ix-1, iz-1)] = W1[Id(iy, ix, iz)] / r[IC(iy-1, ix-1, iz-1)];
	v[IC(iy-1, ix-1, iz-1)] = W2[Id(iy, ix, iz)] / r[IC(iy-1, ix-1, iz-1)];
	w[IC(iy-1, ix-1, iz-1)] = W3[Id(iy, ix, iz)] / r[IC(iy-1, ix-1, iz-1)];
	
	__syncthreads();
	
	k[IC(iy-1, ix-1, iz-1)] = u[IC(iy-1, ix-1, iz-1)]*u[IC(iy-1, ix-1, iz-1)] +
									  v[IC(iy-1, ix-1, iz-1)]*v[IC(iy-1, ix-1, iz-1)] +
									  w[IC(iy-1, ix-1, iz-1)]*w[IC(iy-1, ix-1, iz-1)] ;

	}		
	__syncthreads();	
}
*/


__global__ void KE_kernel(ptype *KE, ptype *KE_reduced)
{
	unsigned int tid   = threadIdx.x;
	unsigned int start = 2*blockIdx.x*blockDim.x;
	__shared__ ptype tile[2*Bsz3D];
	
   if ( start + tid < d_Nc )
        tile[tid] = KE[start + tid];
   else
    	tile[tid] = 0;
	
  
   if ( start + blockDim.x + tid < d_Nc )
        tile[blockDim.x + tid] = KE[start + blockDim.x + tid];
   else
        tile[blockDim.x + tid] = 0;
	
	for(int i = Bsz3D ; i > 0 ; i /= 2)
	{
		__syncthreads();

		if(tid < i)
			tile[tid] += tile[tid + i];
	}
	__syncthreads();		
	
	KE_reduced[blockIdx.x] = tile[0]; 
			
}

__global__ void PfluctuationSqr_Kernel(ptype *d_PflucSqr, ptype *d_p, ptype *d_Pmean)
{
	int ix = blockIdx.x*blockDim.x + threadIdx.x;
	
	if(ix > d_Nc)	{return;}
	ptype Pmean = *d_Pmean;
	
	d_PflucSqr[ix] = (d_p[ix] - Pmean)*(d_p[ix] - Pmean);
	
	__syncthreads();
	
}

/**************************************************************************************************************************
 * -----------------------------------------------------------------------------------------------------------------------*
 **************************************************************************************************************************/
