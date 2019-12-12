function [x,success] = LDLTSolve(A,b)
%ModifiedCholesky Performs a modified Cholesky LDL' Decomposition and
%solves the system
    [n,~]=size(A);
    L = zeros(n,n);
    success=1;
    for i = 1:n
        L(i,i)=1;
    end
    D = zeros(n,1);
    for i=1:n
        D(i)=A(i,i)-L(i,1:i-1).^2*D(1:i-1);
        if isnan(D(i))
            warning('NaN Diagonal element')
            success=0; 
            break
        end
        for j=i+1:n
            L(j,i)=(A(j,i)-L(j,1:i-1).*L(i,1:i-1)*D(1:i-1))/D(i);
            if isnan(L(j,i))
                warning('NaN L element')
                success=0;
                break
            end
        end
        if success==0
            break
        end
    end
if success==1
    y=zeros(n,1);
        y(1)=b(1);
        for i=2:n
            y(i)=b(i)-L(i,1:i-1)*y(1:i-1);
        end
    z=zeros(n,1);
        for i=1:n
            z(i)=y(i)/D(i);
        end
    x=zeros(n,1);
        x(n)=z(n);
        for i=n-1:-1:1
            x(i)=z(i)-L(i+1:n,i)'*x(i+1:n);
        end
else
    x=zeros(size(D));

end


