% the phase field equation with Neumann boundary cindition.
% return  ~~~~~~~~~~~~~~~~~~~~~~~~
% eigvTs_phase  total eigenvalue of the phase field 
% tempux_phase  eigenfunctions in x-direction
% tempuy_phase  eigenfunctions in y-direction
if (MHD_bc=='N')
    [eigfx_phase,eigvx_phase]=eig(Dxx,bmxc,'chol');             % x-direction
    [lambdax_phase,I]=sort(diag(eigvx_phase));
    lambdax_phase(1)=0;
    eigfx_phase=eigfx_phase(:,I);
    [eigfx_phase,~]=Gram_Schmidt(eigfx_phase,diag(bmxc));

    [eigfy_phase,eigvy_phase]=eig(Dyy,bmyc,'chol');             % y-direction
    [lambday_phase,I]=sort(diag(eigvy_phase));
    lambday_phase(1)=0;
    eigfy_phase=eigfy_phase(:,I);
    [eigfy_phase,~]=Gram_Schmidt(eigfy_phase,diag(bmyc));
elseif(MHD_bc=='P')
    Ipx = eye(size(bmxc)); Ipx(1,end) = 1; Ipx(end,1) = 1;
    bmxcp = Ipx'*bmxc*Ipx; Dxxp = Ipx'*Dxx*Ipx;
    bmxcds = cifx *bmx*cifx';
    %%
    Dxxp2 = Dxxp(1:end-1,1:end-1); bmxcp2 = bmxcp(1:end-1,1:end-1);
    [eigfx_phase,eigvx_phase]=eig(Dxxp2,bmxcp2,'chol');             % x-direction
    [lambdax_phase,I]=sort(diag(eigvx_phase));
    lambdax_phase(1)=0;
    eigfx_phase=eigfx_phase(:,I);
    eigfx_phase(end+1,:)=eigfx_phase(1,:);
    [eigfx_phase,~]=Gram_Schmidt(eigfx_phase,diag(bmxc));
    %%
    Ipy = eye(size(bmyc)); Ipy(1,end) = 1; Ipy(end,1) = 1;
    bmycp = Ipy'*bmyc*Ipy; Dyyp = Ipy'*Dyy*Ipy;
    bmycds = cify*bmy*cify';
    %%
    Dyyp2 = Dyyp(1:end-1,1:end-1); bmycp2 = bmycp(1:end-1,1:end-1);
    [eigfy_phase,eigvy_phase]=eig(Dyyp2,bmycp2,'chol');             % x-direction
    [lambday_phase,I]=sort(diag(eigvy_phase));
    lambday_phase(1)=0;
    eigfy_phase=eigfy_phase(:,I);
    eigfy_phase(end+1,:)=eigfy_phase(1,:);
    [eigfy_phase,~]=Gram_Schmidt(eigfy_phase,diag(bmyc));
end


[temp_eigvx,temp_eigvy]=meshgrid(lambdax_phase,lambday_phase);


eigvTs_phase= temp_eigvx'+temp_eigvy';   % the totally eigenvalues 

%--------------------------------------------------------------------~
        [rl,cl]=size(eigvTs_phase);                                % | 
                                                                   % |       
        I1=sparse(eye(rl));I2=sparse(diag(1:cl));                  % |
        transm=sparse(kron(I2,I1));                                % |
               idxeig=diag(transm);                                % |
       tempuy_phase=eigfy_phase(:,idxeig);                         % |         % |
       tempux_phase=repmat(eigfx_phase,1,cl);                      % |         % |
%--------------------------------------------------------------------~
    [nr,nc]=size(eigvTs_phase);
    m_Nx = floor(nr^(betasvv/2)); m_Ny = floor(nc^(betasvv/2));
 %    m_Nx = floor(2*nr/4); m_Ny = floor(nc*2/4);
    Q_x=[zeros(1,m_Nx),...
        exp(-((nr-((m_Nx+1):nr))./(m_Nx-((m_Nx+1):nr)).^2))];
    Q_y=[zeros(1,m_Ny),...
        exp(-((nc-((m_Ny+1):nc))./(m_Ny-((m_Ny+1):nc)).^2))];
    [Qx,Qy] = meshgrid(Q_x,Q_y);  Qxy = (Qx'+Qy')/2;
    eigvTsvv_phase=temp_eigvx'.*Qx'+temp_eigvy'.*Qy';
% test for 

%ut=cos(2*x).*cos(3*y);

%umn=  eigfx_phase'*bmxc*ut*bmyc*eigfy_phase;

%ut2=tempux_phase*sparse(diag(sparse(umn(:))))*tempuy_phase';

