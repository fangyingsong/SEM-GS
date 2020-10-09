function [xmc,xm,D2,D1,bmx,bm,Cinterf,Cinterfdssum,L,DL,xmc2,D2Leg]=genmeo(xg,xi,wi,Di,N)
%inputs: xg the elemets points 
%        xi the GLL points in [-1,1]
%        wi the weight in [-1,1]
%        Di the 1st order derivative in [-1,1]
%        N  the polynomial degree in each element
El = length(xg)-1; %compute the number of element in 1D.
xm = []; % return the mesh grid
xm2 =[]; % return the mesh in 2nd grid
bm = []; % return the mass matrix
D  = []; % return the 1st order derivative matrix
L  = []; % return the interpolation matrix of 1 to 2 
DL = []; % return the 1st order derivative matrix on grid 2
Leg  = []; % for svv matries on Legendre basic
DLeg = [];
Q  = [];   % smooth fliter~ Q^{1/2}!!!!
% solve the pressure on grid 2 we use PN_PN_2 method

N1 = N+1;
%N2 = N1-2;
N2=129;
Cinterf = zeros(El*(N+1)-(El-1),El*(N+1));
Cinterf(1:(N+1),1:(N+1)) = eye(N+1); % the interface  connect information~
Cinterf2= zeros(El*(N2)-(El-1),El*(N2));
Cinterf2(1:(N2),1:(N2)) = eye(N2);

[x2,~]=legslb(N2);        % the pressure grid x2.
L0  = FDMX(0,x2,xi);
DL0 = FDMX(1,x2,xi);
%                 svv in Legendre polynomial
     for j = 0:N
       for i = 0:N
        Leg0(i+1,j+1) = Legendre(j,0,xi(i+1));    %pressure basis
        DLeg0(i+1,j+1)= Legendre(j,1,xi(i+1));
       end
     end   
     mN=floor(sqrt(N)); %   m_N=N^{1/2}; for 
  %    mN=N-3; %   m_N=N^{1/2}; for 
  for i = 1:El
    h = (xg(i+1)-xg(i))/2;
    b = (xg(i+1)+xg(i))/2;
    xm = [xm;xi*h+b];
    xm2= [xm2;x2*h+b];
    bm = blkdiag(bm,sparse(diag(wi*h)));
    D  = blkdiag(D,Di/h);
    %  ~~~~~~~ grid 2 Legendre polynomial
    L   = blkdiag(L,L0);
    DL  = blkdiag(DL,DL0/h);
    %   the pressure basic
    %------------------------------------------this add for svv
    Leg  = blkdiag(Leg,Leg0);
    DLeg = blkdiag(DLeg,DLeg0/h);
        
     for k=0:N
         if k<=mN
            q(k+1)=0;
         else
            q(k+1)=sqrt(exp(-((N-k)/(mN-k))^2));
         end
     end
     Q0=diag(q);
     Q=blkdiag(Q,Q0);
   % clear Leg0 DLeg0
  end
    clear L0 DL0   Leg0 DLeg0 Q0 
 %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  for i = 2:El % the interface dssum
      in1 = (i-1)*(N1)+1;
      in2 = (i-1)*(N1)-(i-1)+1;
      Cinterf(in2:(in2+N),in1:(in1+N)) = eye(N1);
      in1 = (i-1)*(N2)+1;
      in2 = (i-1)*(N2)-(i-1)+1;
      Cinterf2(in2:(in2+N2-1),in1:(in1+N2-1)) = eye(N2);
  end
  Cinterfdssum = diag(1./sum(Cinterf'))*Cinterf;
  Cinterf2dssum = diag(1./sum(Cinterf2'))*Cinterf2;
  D2 = Cinterf*D'*bm*D*Cinterf';  % 2nd order derivative in weak form.
  xmc = Cinterfdssum*xm;
  xmc2= Cinterf2dssum*xm2;
  bmx = Cinterf*bm*Cinterf';
  D1  = D*Cinterf';   % 1st order derivative
  %
  L  = L*Cinterf2';
  % the interpolation form first grid to second grid~
  DL = DL*Cinterf2';
%  ~~~~~~~  Svv matrix----------------------------------------song add
  D1Leg= DLeg*Q*(Leg'*bm*Leg\Leg'*bm)*Cinterf';
  D2Leg= D1Leg'*bm*D1Leg;
end