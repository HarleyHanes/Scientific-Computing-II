function [eigA,iter] = QRalgorithm(A,itermax,tol)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
[Q,R]=qr(A); S=Q;
iter=1;
while 1
    Aold=A;
    A=R*Q;
    [Q,R]=qr(A);
    S=S*Q;
%     err=0;
%     for i=1:length(A)-1
%         err=err+abs(A(i,i+1:end)); %Check abs of lower triag elements
%     end
    err=norm(tril(A,-1))/norm(A); %calculate error as relative power of lower triangular section
    %err=(A*ones(size(A,1),1)-Aold*ones(size(A,1),1));
        if err <= tol
            fprintf("QR Algorithm Successfully converged\n")
            success=1;
            break
        elseif iter==itermax
            fprintf("QR Algorithm Unsuccessful, err=%.2e\n",err)
            success=0;
            break
        end
    iter=iter+1;
end
eigA=diag(A);
end

