%general the meshes of 2D        ------------------------------------------
lagdt=0; 
% MHD_bc magnetic equation boundary conditions.
VELx_bc='P'; VELy_bc='P';    
% velocity boundary conditions 
% 'P' == periodic;  'V' == Dirichlet ; 'N' == Neumann;
for dt=1e-1
alphasvv=1.0;
[x,y]=meshgrid(xmc,ymc);
x=x';y=y'; 

%--------------------------------------------------------------------------

vel = cell(2,1);

vel{1} = ones(size(x));
vel{2} = zeros(size(x));

[m,n]=size(x);  % index for the positions

    for i=1:m
    for j=1:n
        x0=x(i,j);
        y0=y(i,j);
%        if(x0>-0.25&x0<0.25&y0>0.25&y0<0.75)
         if(abs(x0)<=0.08&&abs(y0)<0.08)
            vel{1}(i,j)=1/2;
            vel{2}(i,j)=1/4;
         end
    end
    end
    
    veln = vel;         % u^{n}
    veln1 = vel;        % u^{n-1}
    
    velns{1} = (3*veln{1}-veln1{1})./2;        % u^{*,n+1/2}
    velns{2} = (3*veln{2}-veln1{2})./2;        % u^{*,n+1/2}
% eigen decompsition for velocity equation

eigen_velocity;      % EVP for velocity eq
% return back
% eigvTs_vel  total eigenvalue of the velocity field 
% tempux_vel  eigenfunctions in x-direction
% tempuy_vel  eigenfunctions in y-direction
%----------------initialize the coefficients
T=10000;  idt=1/dt;  simt=0;  % time and time step 
Nstep=floor(T/dt);   % total time steps
ru = 2e-5; rv = ru/2; 
% type alpha
% c = 0.01; kappa = 0.047;
alpha = 1.0;
typeind=1;
if (typeind==1)
% type rho1
c = 0.09; kappa = 0.059; typ='rho1';
elseif(typeind==2)
% type rho2
c = 0.102; kappa = 0.055; typ ='rho2';
elseif(typeind==3)
% type sigma
c = 0.09; kappa = 0.057; typ = 'sigma1';
elseif(typeind==4)
% type sigma2
c = 0.11; kappa = 0.0523; typ = 'sigma2';
elseif(typeind==5)
% type theta
c = 0.03; kappa = 0.057; typ = 'theta1';
elseif(typeind==6)
% type theta2
c = 0.038; kappa = 0.061; typ = 'theta2';
elseif(typeind==7)
  c = 0.014; kappa = 0.053; typ = 'alpha2';
elseif(typeind==8)
% type theta2
c = 0.022; kappa = 0.051; typ = 'gamma1';
elseif(typeind==9)
  c = 0.026; kappa = 0.055; typ = 'gamma2';
elseif(typind==10)
  c = 0.01; kappa = 0.047; typ = 'alpha1';
end
% temperal variables 
v14=zeros(size(x));

filep=strcat('GSM_',typ,'/');
%         save(strcat(filep,'x','.txt'),'x','-ascii');
%         save(strcat(filep,'y','.txt'),'y','-ascii');
        
for nst=1:Nstep
    simt=simt+dt; 
    
    %% first step~~~(t_n,t_{n+1/4}];
    v14 = veln{2}+velns{1}.*velns{2}.^2*dt/2;
    u14 = veln{1}+veln{2}-v14;
    %% second step~~(t_{n+1/4},t_{n+3/4}];
    %  start to solve fractional diffusion equations:
    cmn_vel{1}=eigvTs_vel.^(alpha).*(eigfx_vel'*bmxc*u14*bmyc*eigfy_vel);
    cmn_vel{2}=eigvTs_vel.^(alpha).*(eigfx_vel'*bmxc*v14*bmyc*eigfy_vel);
    
    FXn = idt*u14+c.*(1-vel{1}./2);
    FYn = idt*v14-(c+kappa).*vel{2}./2;
    
    cmnx = eigfx_vel'*bmxc*FXn*bmyc*eigfy_vel-ru/2*cmn_vel{1};
    cmny = eigfx_vel'*bmxc*FYn*bmyc*eigfy_vel-rv/2*cmn_vel{2};
    
    uc= cmnx./(idt + c./2 + ru/2*eigvTs_vel.^(alpha));
    vc= cmny./(idt + (c+kappa)./2 + rv/2*eigvTs_vel.^(alpha));
    
    % convert to the physic space 
    %u34=tempux_vel*sparse(diag(sparse(uc(:))))*tempuy_vel';
    %v34=tempux_vel*sparse(diag(sparse(vc(:))))*tempuy_vel';
    u34 = spectal2physical(uc,tempux_vel,tempuy_vel);
    v34 = spectal2physical(vc,tempux_vel,tempuy_vel);
    
    %% third step 
    vel{2} = v34+velns{1}.*velns{2}.^2*dt/2;
    vel{1} = u34+v34-vel{2};
    
    %% end 
    %% update u v to the next time step!!!
    
    veln1 = veln;        % u^{n-1}
    veln  = vel;         % u^{n}
    
    velns{1} = (3*veln{1}-veln1{1})./2;        % u^{*,n+1/2}
    velns{2} = (3*veln{2}-veln1{2})./2;        % u^{*,n+1/2}
    
    if (mod(nst,200)==0||nst==1)
        L8v = max(max(abs(vel{1})))-min(min(abs(vel{1})));
        fprintf('typeind=%d, step=%d, L8-norm=%e\n',typeind,nst,L8v);
        velx=vel{1};
        if ~exist(filep,'dir')
            mkdir(filep);
        end
%         fig=strcat(filep,'velu',int2str(nst),'.txt');
%         save(fig,'velx','-ascii');
%         vely=vel{2};
%         fig=strcat(filep,'velv',int2str(nst),'.txt');
%         save(fig,'vely','-ascii');
         lagdata.vel = vel;
         lagdata.x = x;
         lagdata.y =y;
         fig=strcat(filep,'lagdata',int2str(nst),'.mat');
         save(fig,'lagdata');
    end
    
end

end
