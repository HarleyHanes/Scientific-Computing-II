function [Data] = GenerateHalton(f,dim,numpoints)
%GenerateHalton Generates Halton spaced x points and their funciton values
%   Utilizes MatLab's haltonsetfunction

%Generate x points
p=haltonset(2);
x=net(p,numpoints);
Data.x=(x).*(dim(2,:)-dim(1,:))+dim(1,:);

%Evaluate y data
Data.y=f(Data.x(:,1),Data.x(:,2));

end

