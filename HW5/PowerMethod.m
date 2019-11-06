function [lambda,q] = PowerMethod(A,q0,itermax,tol)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
iter=1; q=q0/norm(q0); lambda=q'*A*q;
    while 1
         q=A*q/norm(A*q);
         lambda(iter+1)=q'*A*q;
         err=norm(lambda(iter+1)-lambda(iter));
        %Set exit conditions
        if err <= tol
            fprintf("Power Method Successfully converged to lambda=%.2f%+.2fi\n",lambda(end),imag(lambda(end)))
            break
        elseif iter==itermax
            fprintf("Power Method Unsuccessful, err=%.2e\n",err)
            break
        end
        iter=iter+1;
    end
end