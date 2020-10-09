function sol=spectal2physical(coeff,eigx,eigy)
    
    coeff=coeff(:);
    cl=size(eigx,1);
    coeff=repmat(coeff',cl,1);
    sol=(eigx.*coeff)*eigy';
end