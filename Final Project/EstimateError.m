function [Err] = EstimateError(numpoints,dim,f,RBF,metric)
%EstimateError Solves for a number of error metrics for the synthetic
%simulations
x=linspace(dim(1,1),dim(2,1),numpoints(1));
y=linspace(dim(1,2),dim(2,2),numpoints(2));
fMesh=NaN(numpoints);
RBFMesh=NaN(numpoints);
for i=1:numpoints(1)
    for j=1:numpoints(2)
        fMesh(i,j)=f(x(i),y(j));
        RBFMesh(i,j)=RBF([x(i) y(j)]);
    end
end
switch metric
    case 'Max'
        Err=max(max(abs(fMesh-RBFMesh)));
    case 'Mean'
        Err=mean(mean(abs(fMesh-RBFMesh)));
    case 'all'
        Err=reshape(abs(fMesh-RBFMesh),1,numpoints(1)*numpoints(2));
end
end

