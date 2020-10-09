% Differentiation matrices for polynomial basis functions;
% alpha: the corresponding nodes of basis function;
% node: need to be calculated in differentiation matrices;
% the order of differentiation;

function Dm = FDMX(M,alpha,node)
N = length(alpha)-1;
K = length(node)-1;
c = zeros(M+1,N+1,N+1);
for k = 0:K
    c(1,1,1) = 1;c1 = 1;
    for n = 1:N
        c2 = 1;
        for v = 0:n-1
            c3 = alpha(n+1)-alpha(v+1);c2 = c2*c3;
            for m = 0:M
                if m == 0
                    c(m+1,n+1,v+1) = (alpha(n+1)-node(k+1))*c(m+1,n,v+1)/c3;
                else
                    c(m+1,n+1,v+1)=((alpha(n+1)-node(k+1))*c(m+1,n,v+1)-m*c(m,n,v+1))/c3;
                end
            end
        end
        for m = 0:M
            if m == 0
                c(m+1,n+1,n+1) = c1*(-alpha(n)+node(k+1))*c(m+1,n,n)/c2;
            else
                c(m+1,n+1,n+1) = c1*(c(m,n,n)-(alpha(n)-node(k+1))*c(m+1,n,n))/c2;
            end
        end
        c1 = c2;
    end
    for j = 0:N
        Dm(k+1,j+1) = c(M+1,N+1,j+1);
    end
end