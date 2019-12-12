function [beta,cond,tcost,success] = SolveRBFPlanarInterp(Data,basisFunc,reg)
%SolveRBFPlanarInterp Outputs the coeffecients, condition number and
%computation time for deriving coeffecients given basis functions,
%reference points and data for RBF Interpolant

%Check Dimensions 
%Data
[npoints,dim]=size(Data.x);

if length(Data.y)~=npoints
    warning("Number of samples in x and y don't match")
end
%Define x
xi=Data.x';
%Allocate A
A=NaN(npoints,npoints+dim+1);
for i=1:length(Data.x)
    A(i,:)=[basisFunc(Data.x(i,:),xi) Data.x(i,:) 1];
end
A(end+1:end+dim+1,:)=[[Data.x'; ones(1,length(Data.x))] zeros(dim+1)];

if reg==1
    A=A+eye(size(A))*10^(-8);
end

%Solve for coefficients
y=[Data.y;zeros(dim+1,1)];
tic
[beta,success]=LDLTSolve(A,y);
%beta=A\y;
tcost=toc;
%Check cond
cond=1/rcond(A'*A);


end

