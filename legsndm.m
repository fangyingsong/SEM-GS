 function x=legsndm(n)

%  x=legs(n) returns n Legendre-Gauss points
%  Eigenmethod is used for computing nodes. 
%  See Formula (3.184) of the book: J. Shen, T. Tang and L. Wang, Spectral Methods:
%  Algorithms, Analysis and Applications, Springer Series in Compuational
%  Mathematics, 41, Springer, 2011.     
%  Last modified on August 30, 2011

 j=[1:n-1];                  % indices    
 A=diag(j./sqrt((2*j-1).*(2*j+1)),1)...
  +diag(j./sqrt((2*j-1).*(2*j+1)),-1);   %  Create Jacobi matrix
 x= sort(eig(sparse(A)));                %  Compute eigenvalues






