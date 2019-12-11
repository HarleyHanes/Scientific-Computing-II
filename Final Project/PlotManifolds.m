%PlotManifolds
%Creates 3D plot of intersecting manifolds to illustrate why you add
%reference points 
%4=(x-1)^2+y^2+z^2
f=@(x,y)(x.^2<=1).*(y.^2<=1).*sqrt(4-(x-1).^2-y.^2);
x=linspace(-2,2,100);
y=linspace(-2,2,100)';
Z=f(x,y);
mesh(x,y,f(x,y))
%4=(x+1)^2+y^2+z^2
g=@(x,y)sqrt(4-(x+1).^2+y.^2);