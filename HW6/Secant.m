function [x] = Secant(x0,x1,f,itermax)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
x(1)=x0;
x(2)=x1;
for i=3:itermax
    x(i)=x(i-1)-f(x(i-1))*(x(i-1)-x(i-2))/(f(x(i-1))-f(x(i-2)));
    if isnan(x(end))
        x(end)=[];
        break
    end
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

