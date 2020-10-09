clear;clc;
close all
%addpath ../
format long e
%xg=[-2,-1,-0.75,-0.5,-0.25,0,0.25,0.5,1:1:9,10];
h=2;
xg=-1:h:1;
%yg=-1:1/2:1;
yg=-1:h:1;
lagN = 0;
for N=[128,10,20,30,40,50,60,70,80,90,100,110]
%for N = 128
% show the geometry ~~~~~~~|
% show_geo;  %~~~~~~~~~~~~~~~| for test
% end show ~~~~~~~~~~~~~~~~| 
[xi,wi] = legslb(N+1); % the GLL points and weights~[-1,1]
D=FDMX(1,xi,xi); % derivative in 1D
% generate the mesh grid and matrices 
[xmc,xm,Dxx,D1x,bmxc,bmx,cifx,cifdx,Lx,DLx,xmc2,D2Legx]...
    = genmeo(xg,xi,wi,D,N); %x-direction
% xmc grid of each subdomain 
% xm grid of the global domain 
% Dxx the matrix of the -Laplace
% D1x the matrix of dx
% mass matries bmxc bmx
% cifdx the interface dssum operator
% Lx the interpolation operator in 2nd grid 
% DLx the 1st order derivative matrix in 2nd grid
[ymc,ym,Dyy,D1y,bmyc,bmy,cify,cifdy,Ly,DLy,ymc2,D2Legy]...
    = genmeo(yg,xi,wi,D,N);  % y-direction

% solve the imcompressible N-S eq.
%SOLVE_ACMEQ;
SOLVE_GSMEQ_accurate;
lagN=lagN+1;
if (lagN==1)
   ulag = vel{1};
   vlag = vel{2};
else
   L2error(lagN-1,1)=N;
   udiff = cifdx*Lx*ulag*Ly'*cifdy'-vel{1};
   L2error(lagN-1,2)=sqrt(sum(sum(bmxc*(udiff).^2*bmyc)));
   [L2,~,H2] = H2norm(udiff,eigfx_vel,bmxc,bmyc,eigfy_vel,eigvTs_vel);
   L2error(lagN-1,3)=H2;
   fprintf('L2dir=%e, L2eig=%e\n',L2error(lagN-1,2), L2);
   vdiff = cifdx*Lx*vlag*Ly'*cifdy'-vel{2};
   L2error(lagN-1,4)=sqrt(sum(sum(bmxc*(vdiff).^2*bmyc)));
   [L2,~,H2] = H2norm(vdiff,eigfx_vel,bmxc,bmyc,eigfy_vel,eigvTs_vel);
   L2error(lagN-1,5)=H2;
    fprintf('N=%d, lag-N=%d\n',N, lagN);
save('L2errorH2dtem3','L2error')
end
end

% cmnx = eigfx_vel'*bmxc*FXn*bmyc*eigfy_vel
function [L2,H1,H2] = H2norm(u,eigfx_vel,bmxc,bmyc,eigfy_vel,eigvTs_vel)
       umn = eigfx_vel'*bmxc*u*bmyc*eigfy_vel;
       H2 = sqrt(sum(sum(umn.^2.*eigvTs_vel.^2)));
       H1 = sqrt(sum(sum(umn.^2.*eigvTs_vel)));
       L2 =  sqrt(sum(sum(umn.^2)));
end
% end solve~