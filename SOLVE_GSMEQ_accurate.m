%general the meshes of 2D        ------------------------------------------
lagdt=0;   %
% MHD_bc magnetic equation boundary conditions.
VELx_bc='P'; VELy_bc='P';    
% velocity boundary conditions 
% 'P' == periodic;  'V' == Dirichlet ; 'N' == Neumann;
%for dt=[1e-3,1/4,1/8,1/16,1/32,1/64,1/128,1/256,1/512]
for dt = 1e-3
    if lagN==0
        dt = 1e-3;
    else
        dt = 1/512;
    end
alphasvv=1.0;
[x,y]=meshgrid(xmc,ymc);
x=x';y=y'; 

%--------------------------------------------------------------------------

vel = cell(2,1);

vel{1} = sin(pi*x).*cos(2*pi*y);
vel{2} = cos(pi*x).*sin(2*pi*y);

[m,n]=size(x);  % index for the positions

%     for i=1:m
%     for j=1:n
%         x0=x(i,j);
%         y0=y(i,j);
% %        if(x0>-0.25&x0<0.25&y0>0.25&y0<0.75)
%          if(abs(x0)<=0.06&&abs(y0)<0.06)
%             vel{1}(i,j)=1/2;
%             vel{2}(i,j)=1/4;
%          end
%     end
%     end
    
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
T=1.0;  idt=1/dt;  simt=0;  % time and time step 
Nstep=floor(T/dt);   % total time steps
ru = .005; rv = .005/2; c = 0.03; kappa = 0.055;
alpha = 1.0;
% temperal variables 
v14=zeros(size(x));

filep='GSM_Fk61/';
%         save(strcat(filep,'x','.txt'),'x','-ascii');
%         save(strcat(filep,'y','.txt'),'y','-ascii');
        
for nst=1:Nstep
    simt=simt+dt; 
    
    %% first step~~~(t_n,t_{n+1/4}];
%     v14 = veln{2}+velns{1}.*velns{2}.^2*dt/2;
%     u14 = veln{1}+veln{2}-v14;
    %% second step~~(t_{n+1/4},t_{n+3/4}];
    %  start to solve fractional diffusion equations:
    cmn_vel{1}=eigvTs_vel.^(alpha).*(eigfx_vel'*bmxc*vel{1}*bmyc*eigfy_vel);
    cmn_vel{2}=eigvTs_vel.^(alpha).*(eigfx_vel'*bmxc*vel{2}*bmyc*eigfy_vel);
    
    FXn = idt*vel{1}+c.*(1-vel{1}./2)-velns{1}.*velns{2}.^2;
    FYn = idt*vel{2}-(c+kappa).*vel{2}./2+velns{1}.*velns{2}.^2;
    
    cmnx = eigfx_vel'*bmxc*FXn*bmyc*eigfy_vel-ru/2*cmn_vel{1};
    cmny = eigfx_vel'*bmxc*FYn*bmyc*eigfy_vel-rv/2*cmn_vel{2};
    
    uc= cmnx./(idt + c./2 + ru/2*eigvTs_vel.^(alpha));
    vc= cmny./(idt + (c+kappa)./2 + rv/2*eigvTs_vel.^(alpha));
    
    % convert to the physic space 
    %u34=tempux_vel*sparse(diag(sparse(uc(:))))*tempuy_vel';
    %v34=tempux_vel*sparse(diag(sparse(vc(:))))*tempuy_vel';
    vel{1} = spectal2physical(uc,tempux_vel,tempuy_vel);
    vel{2} = spectal2physical(vc,tempux_vel,tempuy_vel);
    
%     %% third step 
%     vel{2} = v34+velns{1}.*velns{2}.^2*dt/2;
%     vel{1} = u34+v34-vel{2};
    
    %% end 
    %% update u v to the next time step!!!
    
    veln1 = veln;        % u^{n-1}
    veln  = vel;         % u^{n}
    
    velns{1} = (3*veln{1}-veln1{1})./2;        % u^{*,n+1/2}
    velns{2} = (3*veln{2}-veln1{2})./2;        % u^{*,n+1/2}
    
%     if (mod(nst,200)==0)
%         nst
%         velx=vel{1};
%         if ~exist(filep,'dir')
%             mkdir(filep);
%         end
%         fig=strcat(filep,'velu',int2str(nst),'.txt');
%         save(fig,'velx','-ascii');
%         vely=vel{2};
%         fig=strcat(filep,'velv',int2str(nst),'.txt');
%         save(fig,'vely','-ascii');
%     end
    
end
lagdt = lagdt+1;
% if (lagdt==1)
%    ulag=veln{1};
%    vlag=veln{2};
% else
%    fprintf('dt=%e, lag-time=%d\n',dt, lagdt);
%    L2error(lagdt-1,1) = dt;
%    udiff = ulag-vel{1};
%    L2error(lagdt-1,2) = sqrt(sum(sum(bmxc*(udiff).^2*bmyc)));
%     [~,~,H2] = H2norm(udiff,eigfx_vel,bmxc,bmyc,eigfy_vel,eigvTs_vel);
%    L2error(lagdt-1,3) = H2 ;
%    vdiff = vlag-vel{2};
%    L2error(lagdt-1,4) = sqrt(sum(sum(bmxc*(vdiff).^2*bmyc)));
%     [~,~,H2] = H2norm(vdiff,eigfx_vel,bmxc,bmyc,eigfy_vel,eigvTs_vel);
%    L2error(lagdt-1,5) = H2;
%     save('L2errordtN128','L2error')
% end

end


% function [L2,H1,H2] = H2norm(u,eigfx_vel,bmxc,bmyc,eigfy_vel,eigvTs_vel)
%        umn = eigfx_vel'*bmxc*u*bmyc*eigfy_vel;
%        H2 = sqrt(sum(sum(umn.^2.*eigvTs_vel.^2)));
%        H1 = sqrt(sum(sum(umn.^2.*eigvTs_vel)));
%        L2 =  sqrt(sum(sum(umn.^2)));
% end
