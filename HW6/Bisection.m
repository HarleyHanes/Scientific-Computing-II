function [x] = Bisection(x0,x1,f,tol)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
if f(x0)>0 && f(x1)<0
    a=x0;
    b=x1;
elseif f(x0)<0 && f(x1)>0
    b=x0;
    a=x1;
else
    warning("f(x0)=%.2f, f(x1)=%.2f",f(x0),f(x1))
end
x=[x0 x1]; iter=3;

while 1
    x(iter)=(a+b)/2;
    err=abs(x(iter)-x(iter-1));
    if err < tol
        break
    elseif f(x(iter)) > 0
        a=x(iter);
    elseif f(x(iter)) < 0
        b=x(iter);
    end
    iter=iter+1;
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