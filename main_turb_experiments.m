clear;clc;
close all
addpath ../
format long e
%xg=[-2,-1,-0.75,-0.5,-0.25,0,0.25,0.5,1:1:9,10];
h=2;
xg=-1:h:1;
%yg=-1:1/2:1;
yg=-1:h:1;
lagN = 0;
for N=[312]
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
SOLVE_GSMEQ2;
% SOLVE_FGMEQ2
% lagN=lagN+1;
% if (lagN==1)
%    ulag=vel{1};
% else
%    L2error(lagN-1,1)=N;
%    L2error(lagN-1,2)=sqrt(sum(sum(bmxc*(cifdx*Lx*ulag*Ly'*cifdy'-vel{1}).^2*bmyc)));
% %    L2error(lagdt,3)=sqrt(sum(sum(bmxc*(vel{2}-Lx*vel{1}*Ly').^2*bmyc)));
% save('L2error10','L2error')
% end
end
%  end solve~
