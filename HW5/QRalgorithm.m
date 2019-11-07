function [eigA,err] = QRalgorithm(A,iter)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
[Q,R]=qr(A); S=Q;
for i=1:iter
    A=R*Q;
    [Q,R]=qr(A);
    S=S*Q;
    err=norm(tril(A,-1))/norm(A);
    eigA(:,i)=diag(A);
end

end

