%% Take Home Midterm
%Harley Hanes
%Sections
%-Convergence Rate Analysis: Code for plotting estimation of convergenc
%   rate
%-RCond Analysis: Code for comparing rcond(A) with different exponential
%                   basis and matrix sizes
%-Functions:SymGS- symmetric gauss-seidel

%% Convergence Rate Analysis
RCondition=NaN(50,1); tol=10^(-6); itermax=100; Success=1;
for j=2:75
    n=j;
    [U,R]=qr(rand(n));
    Z=zeros(n);
    b=rand(n,1);
    for i=1:n
        Z(i,i)=2^(-i);
    end
    A=U*Z*U';
    x0=zeros(n,1);
    ConvergenceRate(j)=-log(max(abs(eig(A))));
    MinSingValue(j)=2^(-j);
        if Success==1
            [x,err,iter,Success]=SymGS(A,b,x0,tol,itermax);
        if Success==0
            Failmark=j;
        end
    end
end

loglog(MinSingValue(1:Failmark-1),ConvergenceRate(1:Failmark-1),'b')
hold on
loglog(MinSingValue(Failmark:end),ConvergenceRate(Failmark:end),'r')
xlabel('Minimum Singular Value')
ylabel('Convergence Rate')

legend('Successful Tests','Unsuccessful Tests')
%% RCond Analysis
%Simulation Setup
RCondition=NaN(50,1); tol=10^(-6); itermax=1000000; Success=1;
base=logspace(0,1/3.3,50);

%Variable Size
for j=2:75
    n=j;
    [U,R]=qr(rand(n));
    Z=zeros(n);
    b=rand(n,1);
    for i=1:n
        Z(i,i)=2^(-i);
    end
    A=U*Z*U';
    x0=zeros(n,1);
    if Success==1
        [x,err,iter,Success]=SymGS(A,b,x0,tol,itermax);
        if Success==0
            Failmark=j;
        end
    end
    RCondition(j)=rcond(A);
end

%Plot
subplot(1,2,1)
semilogy(1:Failmark-1,RCondition(1:Failmark-1),'b')
hold on
semilogy(Failmark:75,RCondition(Failmark:n),'r')
xlabel('Matrix Size')
ylabel('RCond Value')

%Variable Base
RCondition=NaN(length(base),1);
for j=1:length(base)
    n=75;
    [U,R]=qr(rand(n));
    Z=zeros(n);
    b=rand(n,1);
    for i=1:n
        Z(i,i)=base(j)^(i);
    end
    A=U*Z*U';
    x0=zeros(n,1);
    if Success==1
        [x,err,iter,Success]=SymGS(A,b,x0,tol,itermax);
        if Success==0
            Failmark=j;
        end
    end
    RCondition(j)=rcond(A);
end
%Plot
subplot(1,2,2)
loglog(base(1:Failmark-1),RCondition(1:Failmark-1),'b')
hold on
loglog(base(Failmark:end),RCondition(Failmark:end),'r')
xlabel('Exponential Base')
legend('Successful Tests','Unsuccessful Tests')





function [x,err,iter,Success]=SymGS(A,b,x0,tol,itermax)
n=size(A);
if length(A)~= length(b) || length(b) ~= length(x0)
    fprintf('Error!! Sizings of A,b,x0 do not match')
    fprintf('Size of A: [%i %i]',size(A))
    fprintf('Length of b: %i', length(b))
    fprintf('Length of x0: %i', length(x0))
end

iter=0; xNew=x0; err=1; errOld=100;
while  err > tol && iter < itermax
    %Forward
    xOld=xNew;
    xForward=NaN(length(x0),1);
    xBackward=NaN(length(x0),1);
    for i=1:n
        LSum=A(i,1:i-1)*xForward(1:i-1);
        USum=A(i,i+1:n)*xOld(i+1:n);
        xForward(i)=1/A(i,i)*(b(i)-LSum-USum);
    end
    %Backward
    for i=n:-1:1           %Do n:1
        LSum=A(i,1:i-1)*xForward(1:i-1);
        USum=A(i,i+1:n)*xBackward(i+1:n);
        xBackward(i)=1/A(i,i)*(b(i)-USum-LSum);
    end
    xNew=xBackward;
    approxerr=norm(xNew-xOld);
    err=norm(approxerr/xNew);
    %%err=norm(A*x-b);
    iter=iter+1;
%     if mod(iter,10000)==0
%         fprintf('Iteration: %i\n',iter);
%         fprintf('Norm of Residual: %e\n',norm(A*xNew-b));
%         fprintf('Norm of Solution Error: %e\n',err);
%     end
end
x=xNew;
    if iter<itermax
        Success=1;
    else
        Success=0;
        fprintf('Solver did not reach tolerance within %i iterations\n',itermax)
    end
end
