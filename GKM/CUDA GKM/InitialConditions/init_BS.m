%%%%%%%Initialization routine for 3D turbulence%%%%%%%
%%%%Use Erlebacher's method%%%%

clear all;close all;

N = 63;
Np = N*N*N;	%Total number of points

kmax = .5*(N-1);
kmin = -kmax;
krange = kmin:kmax;

I = sqrt(-1);		%complex i
%%Parameters for initial conditions%%
rho0  = 1.0;
Mt    = .3;
urms  = 1.;

k0  = 4;
A =  32*sqrt(2/pi)*urms^2/k0^5;%*2
Re = 30;

%%%Generate random data%%%
u = rand(N,N,N);  uks = fft3(u);
v = rand(N,N,N);  vks = fft3(v);
w = rand(N,N,N);  wks = fft3(w);

%%%%Make velocity field divergence free%%%%
kfft = [0:kmax -kmax:-1];
[k1,k2,k3] = meshgrid(kfft,kfft,kfft);
ksq     = max(k1.^2 + k2.^2 + k3.^2,1.e-14);

kdotvel = k1.*uks + k2.*vks + k3.*wks;

uks = uks - kdotvel.*k1./ksq;
vks = vks - kdotvel.*k2./ksq;
wks = wks - kdotvel.*k3./ksq;

%%%% Compute obtained spectrum%%%%
kEmax = 1+ ceil(sqrt(3)*kmax);
Espec  = zeros(1,kEmax);
for k3 = krange 
  for k2 = krange
    for k1 = krange
     %%%%Map k1,k2,k3 to matlab fft coords
 
     k1m = k1 + 1 +N*(k1 < 0);
     k2m = k2 + 1 +N*(k2 < 0);
     k3m = k3 + 1 +N*(k3 < 0);
     k   = sqrt(k1^2 + k2^2 + k3^2); k123 = max(k,1.e-14);
     km   = 1 + round(k); 
     
     Espec(km) = Espec(km) + abs(uks(k1m,k2m,k3m)).^2 + ...
 	                     abs(vks(k1m,k2m,k3m)).^2 + ...
	                     abs(wks(k1m,k2m,k3m)).^2;
    end;
  end;
end;

Espec = max(Espec,1.e-14);
%%%% Scale spectrum to get desired spectrum %%%%
specrange = 0:kEmax-1; 
kk = specrange;
 
E_ex(kk+1) = sqrt(2/pi)*(16/3)*Mt^2*kk.^4.*exp(-2*kk.^2/k0^2)*Np^2/(k0^5);
for k3 = krange
  for k2 = krange
    for k1 = krange
      %%%%Map k1,k2,k3  to matlab fft coords
    
      k1m  = k1 + 1 +N*(k1 < 0);
      k2m  = k2 + 1 +N*(k2 < 0);
      k3m  = k3 + 1 +N*(k3 < 0);
      km   = 1 + round(sqrt(k1^2 + k2^2 + k3^2)); 
    
      fac                 = sqrt(E_ex(km)/Espec(km));
      uks(k1m,k2m,k3m)    = uks(k1m,k2m,k3m)*fac;
      vks(k1m,k2m,k3m)    = vks(k1m,k2m,k3m)*fac;
      wks(k1m,k2m,k3m)    = wks(k1m,k2m,k3m)*fac;
    
    end;
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Do ifft to get original velocities          %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u = ifft3(uks); u = real(u);
v = ifft3(vks); v = real(v);
w = ifft3(wks); w = real(w);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Verify that the obtained spectrum and     %%%
%%%  velocities are consistent with the        %%%
%%%  desired result                            %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[k1,k2,k3] = meshgrid(kfft,kfft,kfft);
dudxhat = I*k1.*uks; dudx_f = ifft3(dudxhat);
dvdyhat = I*k2.*vks; dvdy_f = ifft3(dvdyhat);
dwdzhat = I*k3.*wks; dwdz_f = ifft3(dwdzhat);

i = 1:N; ip = [2:N 1]; im = [N 1:N-1];
j = 1:N; jp = [2:N 1]; jm = [N 1:N-1];
k = 1:N; kp = [2:N 1]; km = [N 1:N-1];

dx = 2*pi/N;
dudx_c = .5*(u(i,jp,k) - u(i,jm,k))/dx;
dvdy_c = .5*(v(ip,j,k) - v(im,j,k))/dx;
dwdz_c = .5*(w(i,j,kp) - w(i,j,km))/dx;

%%%%1)  Check total enery%%%%%
tke = sum( u(:).^2 + v(:).^2 + w(:).^2)/Np; 

uks = uks/Np; vks = vks/Np; wks = wks/Np;
ek  = sum( abs(uks(:)).^2 + abs(vks(:)).^2 + abs(wks(:)).^2); 
      
%%%%2) Check that tke = ek = Integral( E(k)*dk ) %%%
ek_ex = 3*A/64*sqrt(2*pi)*k0^5;
Espec = zeros(size(Espec));
for k3 = krange
  for k2 = krange
    for k1 = krange
     %%%%Map k1,k2 to matlab fft coords
     k   = sqrt(k1^2 + k2^2 + k3^2);
 
     k1m = k1 + 1 +N*(k1 < 0);
     k2m = k2 + 1 +N*(k2 < 0);
     k3m = k3 + 1 +N*(k3 < 0);
     km  = round(k) + 1;
 
     Espec(km) = Espec(km) + abs(uks(k1m,k2m,k3m)).^2 + ...
 	                     abs(vks(k1m,k2m,k3m)).^2 + ...
	  	             abs(wks(k1m,k2m,k3m)).^2;
    end;
  end;
end;
plot(kk,Espec(kk+1),'-o');hold on;
E_ex = E_ex/Np^2;
plot(kk,E_ex(kk+1),'-rx');
  

%%%%%%Compute viscosity coeff %%%%%%%
dvdy      = real(dvdy_f(:));
urms      = sqrt(ek/3);


crms = sqrt(3)*urms/Mt; %sqrt(3)*urms/Mt; %U0 = sqrt(3)*urms;
gam  = 1.4;
p0   = rho0*crms*crms/gam;


lambda1    = urms/sqrt(sum(dvdy.^2)/Np);
lambda = 2./k0;
mu        = urms*lambda*rho0/Re;
eddy_time = lambda/urms;

%%%%%% Write data out %%%%%%
fid = fopen('Init.dat','w');

fwrite(fid, N, 'int');
fwrite(fid, rho0, 'float');
fwrite(fid, Mt   , 'float');
fwrite(fid, Re   , 'float');
fwrite(fid, k0   , 'float');
fwrite(fid, urms , 'float');
fwrite(fid, mu , 'float');
fwrite(fid, crms , 'float');
fwrite(fid, p0 , 'float');
fwrite(fid, eddy_time , 'float');
fwrite(fid, u, 'float');
fwrite(fid, v, 'float');
fwrite(fid, w, 'float');
fclose(fid);
