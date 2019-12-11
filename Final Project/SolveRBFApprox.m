function [beta,cond] = SolveRBFApprox(Data,xi,basisFunc)
%UNTITLED11 Summary of this function goes here
%   Detailed explanation goes here
%Check Dimensions 
%Data
[npoints,dim]=size(Data.x);
[~,NumBasis]=size(xi);
if size(xi,1)~=dim
    warning('xi has dim %i, x has dim %i',size(xi,1),dim)
end
if length(Data.y)~=npoints
    warning("Number of samples in x and y don't match")
end
%Allocate A
A=NaN(npoints,NumBasis);
for i=1:length(Data.x)
    A(i,:)=basisFunc(Data.x(i,:),xi);
end

%Solve for coefficients
beta=(A'*A)\(A'*Data.y);
%Check sizing
if length(beta)~=NumBasis
    warning("Number of coeffecients doesn't match number of basis functions")
end
cond=1/rcond(A'*A);

end

