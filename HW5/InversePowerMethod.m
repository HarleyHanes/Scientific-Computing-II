function [lambda,q,success] = InversePowerMethod(A,q0,mu,itermax,tol,write)
q=q0/norm(q0);
B=A-mu*eye(size(A));
[L,U]=LUDecomp(B);
%[L,U]=lu(B);
iter=1;
lambda=q'*A*q;
while 1
    q=LowerSolve(L,UpperSolve(U,q));
    %q=UpperSolve(U,LowerSolve(L,q));
    q=linsolve(U,linsolve(L,q));
    q=q/norm(q);
    lambda(iter+1)=q'*A*q;
    err=abs(lambda(iter+1)-lambda(iter));
        if err <= tol
            if write
                fprintf("Inverse Power Method Successfully converged to lambda=%.2f%+.2fi\n",lambda(end),imag(lambda(end)))
            end
            success=1;
            break
        elseif iter==itermax
            if write
                fprintf("Inverse Power Method Unsuccessful, err=%.2e\n",err)
            end
            success=0;
            break
        end
    iter=iter+1;
end
end
