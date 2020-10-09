%The code of Lengendre is to evaluate Ln^(m)(x);
function r = Legendre(n,m,x)
s = zeros(n+1,m+1);
for j = 1:m+1
    if j == 1
        s(1,j) = 1;
        s(2,j) = x;
        for k = 2:n
            s(k+1,j) = ((2*k-1)*x*s(k,j)-(k-1)*s(k-1,j))/k;
        end
    else s(1,j) = 0;
        if j == 2
            s(2,j) = 1;
        else s(2,j) = 0;
        end
        for k = 2:n
            s(k+1,j) = (2*k-1)*s(k,j-1)+s(k-1,j);
        end
    end
end
r = s(n+1,m+1);