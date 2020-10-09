% the velocity equation with periodic BC x-direction boundary 
% wall on y-direction boundary 
% return 
% eigvTs_vel  total eigenvalue of the velocity field 
% tempux_vel  eigenfunctions in x-direction
% tempuy_vel  eigenfunctions in y-direction
if (VELx_bc=='P')
Ipx = eye(size(bmxc));   Ipx(1,end) = 1; Ipx(end,1) = 1;
bmxcp = Ipx'*bmxc*Ipx;   Dxxp = Ipx'*Dxx*Ipx;   %~~~~~~~~~~~
bmxcds =cifx*bmx*cifx';
%% for solving the eigenproblems in 1D case.
%   --------   .----.    .--.    .--------.
     Dxxp2=Dxxp(1:end-1,1:end-1); bmxcp2=bmxcp(1:end-1,1:end-1);
     [eigfx_vel,eigvx_vel]=eig(Dxxp2,bmxcp2,'chol');
     [lambdax_vel,I] = sort(diag(eigvx_vel));
     lambdax_vel(1) = 0;
     eigfx_vel = eigfx_vel(:,I);
     eigfx_vel(end+1,:)=eigfx_vel(1,:);
%      
%     Weig=diag(bmxcds);
      Weig=diag(bmxc);
     [eigfx_vel, ~] = Gram_Schmidt(eigfx_vel, Weig);
elseif(VELx_bc=='V')
Dxxd=Dxx(2:end-1,2:end-1); bmxcd= bmxc(2:end-1,2:end-1);   %left D right D
%Dxxd=Dxx(2:end,2:end); bmxcd= bmxc(2:end,2:end);  % left D right Newmann      
[eigfx_vel,eigvx_vel]=eig(Dxxd,bmxcd,'chol');           
[lambdax_vel,I]=sort(diag(eigvx_vel));
eigfx_vel=eigfx_vel(:,I);
eigfx_vel=[eigfx_vel(1,:)*0;eigfx_vel;eigfx_vel(end,:)*0];
%eigfx_vel=[eigfx_vel(1,:)*0;eigfx_vel];


[eigfx_vel,~]=Gram_Schmidt(eigfx_vel,diag(bmxc));

elseif(VELx_bc=='N')
     [eigfx_vel,eigvx_vel]=eig(Dxx,bmxc,'chol');             % x-direction
    [lambdax_vel,I]=sort(diag(eigvx_vel));
    lambdax_vel(1)=0;
    eigfx_vel=eigfx_vel(:,I);
    [eigfx_vel,~]=Gram_Schmidt(eigfx_vel,diag(bmxc));
end   
%--------------------------------------------------------------------------
if (VELy_bc=='P')
Ipy = eye(size(bmyc));   Ipy(1,end) = 1; Ipy(end,1) = 1;
bmycp = Ipy'*bmyc*Ipy;   Dyyp = Ipy'*Dyy*Ipy;   %~~~~~~~~~~~
bmycds =cify*bmy*cify';
%% for solving the eigenproblems in 1D case.
%   --------   .----.    .--.    .--------.
     Dyyp2=Dyyp(1:end-1,1:end-1); bmycp2=bmycp(1:end-1,1:end-1);
     [eigfy_vel,eigvy_vel]=eig(Dyyp2,bmycp2,'chol');
     [lambday_vel,I] = sort(diag(eigvy_vel));
     lambday_vel(1) = 0;
     eigfy_vel = eigfy_vel(:,I);
     eigfy_vel(end+1,:)=eigfy_vel(1,:);
%      
%     Weig=diag(bmycds);
      Weig=diag(bmyc);
     [eigfy_vel, ~] = Gram_Schmidt(eigfy_vel, Weig);
    
elseif(VELy_bc=='V')
Dyyd=Dyy(2:end-1,2:end-1); bmycd= bmyc(2:end-1,2:end-1);
    
[eigfy_vel,eigvy_vel]=eig(Dyyd,bmycd,'chol');           
[lambday_vel,I]=sort(diag(eigvy_vel));
eigfy_vel=eigfy_vel(:,I);
eigfy_vel=[eigfy_vel(1,:)*0;eigfy_vel;eigfy_vel(end,:)*0];

[eigfy_vel,~]=Gram_Schmidt(eigfy_vel,diag(bmyc));


elseif(VELy_bc=='N')
     [eigfy_vel,eigvy_vel]=eig(Dyy,bmyc,'chol');             % y-direction
    [lambday_vel,I]=sort(diag(eigvy_vel));
    lambday_vel(1)=0;
    eigfy_vel=eigfy_vel(:,I);
    [eigfy_vel,~]=Gram_Schmidt(eigfy_vel,diag(bmyc));
end

%%

[temp_eigvx,temp_eigvy]=meshgrid(lambdax_vel,lambday_vel);

eigvTs_vel= temp_eigvx'+temp_eigvy';   % the totally eigenvalues 

%-------------------------------------------------------------------------|
        [rl,cl]=size(eigvTs_vel);                                       % | 
                                                                        % |       
        I1=sparse(eye(rl));I2=sparse(diag(1:cl));                       % |
        transm=sparse(kron(I2,I1));                                     % |
               idxeig=diag(transm);                                     % |
       tempuy_vel=eigfy_vel(:,idxeig);                                  % |
       tempux_vel=repmat(eigfx_vel,1,cl);                               % |
%--------------------------------------------------------------------------
% svv
    [nr,nc]=size(eigvTs_vel);
    m_Nx = floor(nr^(alphasvv/2)); m_Ny = floor(nc^(alphasvv/2));
%    m_Nx = floor(nr*2/4); m_Ny = floor(nc*2/4);
    Q_x=[zeros(1,m_Nx),...
        exp(-((nr-((m_Nx+1):nr))./(m_Nx-((m_Nx+1):nr)).^2))];
    Q_y=[zeros(1,m_Ny),...
        exp(-((nc-((m_Ny+1):nc))./(m_Ny-((m_Ny+1):nc)).^2))];
    [Qx,Qy] = meshgrid(Q_x,Q_y);  Qxy = (Qx'+Qy')/2;   
    eigvTsvv_vel=temp_eigvx'.*Qx'+temp_eigvy'.*Qy';

% test for 

%ut=cos(2*x).*(pi^2-y.^2);

%umn=  eigfx_vel'*bmxc*ut*bmyc*eigfy_vel;

%ut2=tempux_vel*sparse(diag(sparse(umn(:))))*tempuy_vel';
