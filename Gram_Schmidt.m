function [Q, R] = Gram_Schmidt(A, w)
%w is the weight <A_i,A_j>_w=\delta_{i,j}
    [m, n] = size(A);
    Q = complex(zeros(m, n));
    R = complex(zeros(n, n));
    v = zeros(n, 1);

%     for j = 1:n
%         v = A(:,j);
%         for i = 1:j-1
%             R(i,j) = sum(   v    .* conj( Q(:,i) ) .* w ) / ...
%                      sum( Q(:,i) .* conj( Q(:,i) ) .* w );
%             v = v - R(i,j) * Q(:,i);
%         end
%         R(j,j) = norm(v);
%         Q(:,j) = v / R(j,j);
%     end
%     
  for j = 1:n
      v = A(:,j);
      if (j>1)
          if(j==2)
             R(1,j) = sum(   v    .* conj( Q(:,1) ) .* w ) / ...
                      sum( Q(:,1) .* conj( Q(:,1) ) .* w );
             v = v - R(1,j) * Q(:,1); 
          else
      %R(1:j-1,j)=(Q(:,1:j-1)'*diag(w)*conj(Q(:,1:j-1)))\(conj(Q(:,1:j-1)')*diag(w)*v);
         R(1:j-1,j)=(conj(Q(:,1:j-1)')*diag(w)*v);
         v=v-(sum(diag(R(1:j-1,j))*(Q(:,1:j-1))'))';
          end
      end
       R(j,j) = norm(v);
       Q(:,j) = v / R(j,j);
       Q(:,j) = Q(:,j)/sqrt(Q(:,j)'*diag(w)*Q(:,j));
  end
    
end