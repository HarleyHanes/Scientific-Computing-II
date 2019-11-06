%IcholTestCode
%Code for computational efficacy for PCG with Jacobian, Incomplete
%Cholesky, and no preconditioner
%Authors- Harley Hanes, Elliot Hill
clear;clc;close all;
nvec=linspace(20,40,21).^2;
tol=10^(-6); itermax=200;
for i=1:length(nvec)
    [A,b,x0,~]=CreateTestMatrix(nvec(i));
    [~,iterChol(i),~,timeChol(i),convChol(i)]=HW4PCG(A,b,x0,tol,itermax,'Cholesky');
    [~,iterJac(i),~,timeJac(i),convJac(i)]=HW4PCG(A,b,x0,tol,itermax,'Jacobi');
    [~,iterNull(i),~,timeNull(i),convNull(i)]=HW4PCG(A,b,x0,tol,itermax,'None');
end
subplot(2,1,1)
plot(nvec,timeJac)
ylabel('Computation time (s)')
hold on
plot(nvec,timeChol)
plot(nvec,timeNull)
legend('Jacobi','Incomplete Cholesky','No Preconditioner')

subplot(2,1,2)
plot(nvec,iterJac)
hold on
ylabel('Iterations')
plot(nvec,iterChol)
plot(nvec,iterNull)
xlabel('Matrix Size')

figure
plot(nvec,convJac)
hold on
plot(nvec,convChol)
plot(nvec,convNull)
legend('Jacobi','Incomplete Cholesky','No Precondition')
ylabel('$\frac{\sqrt{\kappa}-1}{\sqrt{\kappa}+1}$',"Interpreter","latex")
xlabel('Matrix Size')


function [A,b,x0, ConvergenceRate] = CreateTestMatrix(n)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% A=zeros(n^2);
% for i=1:n^2
%     for j=1:n^2
%         if i==j 
%             A(i,j)=5;
%         elseif j==i-1 || j==i+1 || j==i-10 || j==i+10
%             A(i,j)=rand;
%         end
%     end
% end
[U,S,~] = svd(randn(n));
s = diag(S);
A = U*diag(s+max(s))*U'; % to make A symmetric, well-contioned
for i=1:length(A)
    A(i,i)=sum(abs(A(:,i)))+1; %Make A diagonally dominant
end
b = randn(n,1);
x0= zeros(n,1);
kappa=cond(A);
ConvergenceRate=(sqrt(kappa)-1)/(sqrt(kappa)+1);
% b=rand(length(A),1);
% x0=zeros(length(A),1);
end
