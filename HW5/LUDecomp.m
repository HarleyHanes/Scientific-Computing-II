function [L,U] = LUDecomp(A)
sizecheck=size(A);
if sizecheck(1)~=sizecheck(2)
    fprintf("Error!! Entered matrix is not square")
    keyboard
end
n=sizecheck(1);
U=zeros(n);
L=zeros(n);

for i=1:n
    for k=i:n
        U(i,k)=A(i,k)-L(i,:)*U(:,k);
    end
    for k=i:n
        if i==k
            L(k,k)=1;
        else
            L(k,i)=(A(k,i)-L(k,:)*U(:,i))/U(i,i);
        end
    end
end
%Check LU Decomp
q=rand(size(A,1),1);
err=norm(A*q-L*U*q);
if ~istriu(U) || ~istriu(L') || err>=1e-12
    warning("LU decomposition unsuccessful, err=%.2e",err)
end
end
