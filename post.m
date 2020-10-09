v1e=sin(simt)*sin(2*pi*y).*sin(pi*x).^2*pi;
v2e=-sin(simt)*sin(2*pi*x).*sin(pi*y).^2*pi;
pe=sin(simt-dt/2)*cos(pi*x).*sin(pi*y);
phie=sin(simt)*cos(pi*x).*cos(pi*y);
Phi1e=sin(simt)*cos(y*pi);
%v1e=v1post;
%v2e=v2post;
%pe=L*pnpost*L';
%phie=phipost;
simu_count0=lagdt;
L2nv1(simu_count0)=(sum(sum((bmxc*(vel{1}-v1e).^2*bmyc))))^(1/2);
L2nv2(simu_count0)=(sum(sum((bmxc*(vel{2}-v2e).^2*bmyc))))^(1/2);
L2np(simu_count0)=(sum(sum((bmxc*(cifdx*Lx*pressuren*Ly'*cifdy'...
    -pe).^2*bmyc))))^(1/2);
L2nphi(simu_count0)=(sum(sum((bmxc*(phie-phi).^2*bmyc))))^(1/2);
L2nPhi(simu_count0)=(sum(sum((bmxc*(Phi1e-Phi{1}).^2*bmyc))))^(1/2);


test_omega1=((cifdx*D1x*Phins{1}).*omegans{1}...
    +(Phins{1}*D1y'*cifdy').*omegans{2});
test_omega2=((cifdx*D1x*Phins{2}).*omegans{1}...
    +(Phins{2}*D1y'*cifdy').*omegans{2});
L2nOmega(simu_count0)=(sum(sum((bmxc*(test_omega1./rhons ...
    -omegaEPS{1}(x,y,simt-dt/2)./rhons).^2*bmyc))))^(1/2);
L2nOmega2(simu_count0)=(sum(sum((bmxc*(test_omega2./rhons ...
    -omegaEPS{2}(x,y,simt-dt/2)./rhons).^2*bmyc))))^(1/2);
% splitting error------------------------------
%L2nspv1(simu_count0)=dt*(sum(sum(((FDX*(v10n-v10n1)*FDY).^2))))^(1/2)/4;
%L2nspv2(simu_count0)=dt*(sum(sum(((FDX*(v20n-v20n1)*FDY).^2))))^(1/2)/4;
%L2nspp(simu_count0)=dt*(sum(sum(((D'*bmx*D*(L*psin*L')*D'*bmy*D).^2))))^(1/2);
%L2nsphi(simu_count0)=dt*(sum(sum(((FDX*(u-u0n1)*FDY).^2))))^(1/2)/4;
dtn(simu_count0)=dt;