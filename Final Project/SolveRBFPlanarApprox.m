function [beta,cond,Tcost,success] = SolveRBFPlanarApprox(Data,xi,basisFunc)
%SolveRBFPlanarApprox Outputs the coeffecients, condition number and
%computation time for deriving coeffecients given basis functions,
%reference points and data
%Data
[npoints,dim]=size(Data.x);
[~,NumRef]=size(xi);
if size(xi,1)~=dim
    warning('xi has dim %i, x has dim %i',size(xi,1),dim)
end
if length(Data.y)~=npoints
    warning("Number of samples in x and y don't match")
end
%Allocate A
A=NaN(npoints,NumRef+dim+1);
for i=1:length(Data.x)
    A(i,:)=[basisFunc(Data.x(i,:),xi) Data.x(i,:) 1];
end
%A(end+1:end+dim+1,:)=[[Data.x'; ones(1,npoints)] zeros(dim+1,dim+1)];

%Solve for coefficients
%y=[Data.y;zeros(dim+1,1)];
tic
beta=(A'*A)\(A'*Data.y);      %Linear LSQ
success=1;
Tcost(1)=toc;


% B=A'*A+10^(-8)*eye(size(A'*A));
% 
% tic
%     [beta,success]=LUSolve(B,A'*Data.y);
% Tcost=toc;



%Check sizing
if length(beta)~=NumRef+dim+1
    warning("Number of coeffecients doesn't match number of basis functions")
end

 cond=1/rcond(A'*A);

end