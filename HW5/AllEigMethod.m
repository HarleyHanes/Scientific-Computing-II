function [eigA,success] = AllEigMethod(A,itermax,tol)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
[Q,R]=qr(A); S=Q;
iter=1;
while 1
    Aold=A;
    A=R*Q;
    [Q,R]=qr(A);
    S=S*Q;
    err=(A*ones(size(A),1)-Aold*ones(size(A),1));
        if err <= tol
            fprintf("Inverse Power Method Successfully converged to lambda=%.2f%+.2fi\n",lambda(end),imag(lambda(end)))
            success=1;
            break
        elseif iter==itermax
            fprintf("Inverse Power Method Unsuccessful, err=%.2e\n",err)
            success=0;
            break
        end
end
eigA=diag(A);
end

