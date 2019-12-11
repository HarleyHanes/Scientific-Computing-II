function [A,err] = HessenbergQR
%example
rand('seed',1)
X=rand(7); DD=diag((1:7)); X1=X'*DD*X;  X2=inv(X)*DD*X 
X=X2;

[n,m]=size(X);
if n~=m
    error("A is not square")
end
%A = X;
A=rand(7)+1i*1e-8*eye(7);
P = eye(n);

for k=1:n-2
    %Compute uk and alpha
    x=A(k+1:n,k);
    [u,alpha]= ELhousegen(x);
    
    %Compute PkA
    A(k+1:n,k:n)=A(k+1:n,k:n)-2*u*(u'*A(k+1:n,k:n));
    %Compute PkAPk
    A(1:n,k+1:n)=A(1:n,k+1:n)-2*(A(1:n,k+1:n)*u)*u';
    P(1:n,k+1:n)=P(1:n,k+1:n)-2*(P(1:n,k+1:n)*u)*u';
end
%Check Hessenberg
err=norm(tril(A,-2));
fprintf("Norm of Lower section of phase I: %.2e\n", err)

err=1;iter=1; A1=A,
while iter<100 & err>1e-5
    G=eye(n);
    [A]=HessenToTri(A) ;
    err=max(abs(diag(A,-1)));
    iter=iter+1;
    plot(diag(A),zeros(n,1),'.','MarkerSize',22),
    title(num2str(iter)),axis([-1  10 -2 2]),pause(0.05)
end
fprintf("Norm of Lower section of phase II: %.2e\n", err)
eig(A)
end
%%%%%%%%%%%%%%%%
% Algorithm 12.13 (Generate a Householder transformation) To given
% x ? Rn the following algorithm computes a = ? and the vector u so that
% (I ? uuT )x = ?e1.
function [u,a]= ELhousegen(x)
    rho = -sign(x(1));
    s = norm(x);
    a = rho*s;
    if s==0
    u(1)= 1; return ;
    end

    u = x; u(1) = u(1)-a;
    u = u/sqrt(2*s*(s-rho*x(1)));
end
%%%%%%%%%%%%%%%%
function [H]=HessenToTri(A)

[N M]=size(A);
GT=eye(N);    

H = A;

for i=1:N-1,  ip1=i+1;
    rho=sqrt(H(i,i)^2+H(ip1,i)^2);
    c(i) = H(i,i)/rho;  s(i) = H(ip1,i)/rho;
    H(i:ip1,i:N) = [[c(i) s(i)];[-s(i) c(i)]]*H(i:ip1,i:N);  
    GT(i:ip1,i:N) = [[c(i) -s(i)];[s(i) c(i)]]*GT(i:ip1,i:N); 
end



H1=H;
for i=1:N-1,  ip1=i+1;
    rho=sqrt(H(i,i)^2+H(ip1,i)^2);
    H(1:ip1,i:ip1) = H(1:ip1,i:ip1) *[[c(i) -s(i)];[s(i) c(i)]];  
end
end
% %%%%%%%%%%%%%%%%