function [x,success] = LDLTSolve(A,b)
%ModifiedCholesky Performs a modified Cholesky LDL' Decomposition
%   Based on code example from https://sites.ualberta.ca/~xzhuang/Math381/Lab5.pdf
%       Added checks for numerical success and failure
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
    z = zeros(n,1);
    y = zeros(n,1);
    x = zeros(n,1);
    % Ly=b
        y(1) = b(1);
        for i=2:n
            y(i) = b(i)-L(i,1:i-1)*y(1:i-1);
        end
    % Dz = y;
        for i = 1:n
            z(i)=y(i)/D(i);
        end
    % L^T x = z;
        x(n) = z(n);
        for i = n-1:-1:1
            x(i) = z(i)-L(i+1:n,i)'*x(i+1:n);
        end
else
    x=zeros(size(D));

end

% n=size(A,2);
% L=zeros(n);D=zeros(n,1);
% for i=1:n
%     Sum1=0;
%     for j=1:i-1
%         Sum1=Sum1+L(i,j)^2*D(j);
%     end
%     D(i)=A(i,i)-Sum1;
%     for j=i+1:n
%         Sum2=0;
%         for k=1:i-1
%             Sum2=Sum2+L(j,k)*L(i,k)*D(k);  
%         end
%         L(j,i)=(A(j,i)-Sum2)/D(i);
%     end
% end
% D=diag(D);
%end

