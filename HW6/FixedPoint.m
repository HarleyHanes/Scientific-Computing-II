function [x] = FixedPoint(x0,phi,itermax)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
x=NaN(1,itermax);
x(1)=x0;
for i=2:itermax
    x(i)=phi(x(i-1));
%     err=abs(x-xold);
%     xold=x;
%     if err < tol
%         fprintf("Newton Successfully converged to x=%.2f\n",x)
%         break
%     elseif iter==itermax
%         fprintf("Newton failed to converge\n")
%     end
end
end

