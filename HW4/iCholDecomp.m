function [r] = iCholDecomp(A)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
[n,~]=size(A);
for k=1:n-1
    r(k,k)=sqrt(A(k,k));
    for j=k+1:n
        r(k,j)=A(k,j)/r(k,k);
    end
    for i=k+1:n
        for j=1:n
            if A(i,j)~=0
                A(i,j)=A(i,j)-r(k,i)*r(k,j);
            end
        end
    end
end
r(n,n)=sqrt(A(n,n));
end

