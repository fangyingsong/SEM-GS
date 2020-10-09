function     gphi=gfun(phi,eps)
       gphi=1/eps^2*(phi.^3-phi);
  % gphi=1/eps^2*(phi.^3/4+3*phi.^2.*(phi/4-1/4));
end