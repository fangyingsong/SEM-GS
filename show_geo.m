%show the geometry domian~
[x,y] = meshgrid(xg,yg);
idx = (x<-pi&y<-pi);
idx2 = find(idx);
x2 = x(:);
y2 = y(:);
x2(idx2) = NaN;
y2(idx2) = NaN;
x = reshape(x2,size(x));
y = reshape(y2,size(x));
%x(idx,idx)=[];
%y(idx,idx)=[];
surf(x,y,y*0)
axis equal
clear x2 y2 x y idx idx2 