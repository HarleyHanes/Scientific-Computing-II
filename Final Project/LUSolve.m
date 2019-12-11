function [x,success] = LUSolve(A,b)
sizecheck=size(A);
if sizecheck(1)~=sizecheck(2)
    fprintf("Error!! Entered matrix is not square")
    keyboard
end
n=sizecheck(1);
U=zeros(n);
L=zeros(n);

%Pivot A
Pnet=eye(n);
Ad=A;
for i=1:n
    A0=A;
    P=eye(n);
    maxVal=max(A(i:end,i));
    for j=i:n
        if A(j,i)==maxVal
            A(i,:)=A(j,:);
            A(j,:)=A0(i,:);
            P(i,i)=0;P(j,j)=0;P(i,j)=1;P(j,i)=1;
            Pnet=P*Pnet;
            break
        end
    end
end

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
if ~istriu(U) || ~istriu(L') || err>=1e-10
    warning("LU decomposition unsuccessful, err=%.2e",err)
    success=0;
else
    success=1;
end
%x=LowerSolve(L,UpperSolve(U,b));
x=UpperSolve(U,LowerSolve(L,b));
x=Pnet^(-1)*x;
end
