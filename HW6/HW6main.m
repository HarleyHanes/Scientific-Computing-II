% HW 6 Main
clear; clc; close all;
%% Problem 1-4
    %Define f
    f=@(x)(x+1).*exp((-x.^2/12))-1;
    fp=@(x)1/6*exp(-x.^2/12).*(-x-3).*(x-2);
%     plot(linspace(-1,8,40),f(linspace(-1,8,40)))
    figure

    x0=4; itermax=20;
    %% Newton
    xNew=Newton(x0,f,fp,itermax);
    alpha=xNew(end);
    plot(log(abs(xNew(1:end-2)-alpha)),log(abs(xNew(2:end-1)-alpha)),"*-")
    hold on 
    %% Double Newton
    xDNew=DoubleNewton(x0,f,fp,itermax);
    alpha=xDNew(end);
    plot(log(abs(xDNew(1:end-2)-alpha)),log(abs(xDNew(2:end-1)-alpha)),"*-")

    %% Secant
    xSec=Secant(x0,x0+2,f,itermax);
    alpha=xSec(end);
    plot(log(abs(xSec(1:end-2)-alpha)),log(abs(xSec(2:end-1)-alpha)),"*-")

    %% Bisection
    xBi=Bisection(x0,x0+2,f,10^(-8));
    alpha=xBi(end);
    plot(log(abs(xBi(1:end-2)-alpha)),log(abs(xBi(2:end-1)-alpha)),"*-")


    xlabel("$log|x_k-x_N|$",'Interpreter','Latex')
    ylabel("$log|x_{k+1}-x_N|$",'Interpreter','Latex')
    plot(.5*linspace(-25,.3,30),linspace(-25,.3,30)-.2,"--")
    plot((1/3)*linspace(-25,.3,30),linspace(-25,.3,30)-.2,"--")
    legend('Newton','Double Newton',"Secant","Bisection",'p=2','p=3','Location','Northwest')
    axis([-25 .5 -25 .5])
    % x=linspace(-3,3,200);
    % plot(x,f(x))
%% Problem 5-6
clear; clc;
x0=8; itermax=50;
g=@(x).2*sin(x)+.5;
[x] = FixedPoint(x0,g,itermax);

clear;N=25;x0=1;
g1=@(x) -3/10*(10*x.^3-15*x.^4+6*x.^5)+13/20;
g2=@(x) -4/10*(10*x.^3-15*x.^4+6*x.^5)+14/20;
g3=@(x) -5/10*(10*x.^3-15*x.^4+6*x.^5)+15/20;
g4=@(x) -6/10*(10*x.^3-15*x.^4+6*x.^5)+16/20;
g1p=@(x) -3/10*(30*x.^2-60*x.^3+30*x.^4);
g2p=@(x) -4/10*(30*x.^2-60*x.^3+30*x.^4);
g3p=@(x) -5/10*(30*x.^2-60*x.^3+30*x.^4);
g4p=@(x) -6/10*(30*x.^2-60*x.^3+30*x.^4);
gp={g1p(1/2);g2p(1/2);g3p(1/2);g4p(1/2)};


x(1,:)=FixedPoint(x0,g1,N);
x(2,:)=FixedPoint(x0,g2,N);
x(3,:)=FixedPoint(x0,g3,N);
x(4,:)=FixedPoint(x0,g4,N);
%plot(x(:,1:end-1)',x(:,2:end)')

for i=1:4
    TableCell{i,1}=sprintf("g%i",i);
    TableCell{i,2}=abs(x(i,end)-.5);
end
TableCell(:,3)={abs(g1p(1/2));abs(g2p(1/2));abs(g3p(1/2));abs(g4p(1/2))};
cell2table(TableCell)

%% Problem 7
clear;
D=@(x)1./(x.^2+1);
S=@(x)x+.2*x.^2;
f=@(x)D(x)-S(x);
x=Secant(0,10,f,50);
%plot(1:length(x),x)
fprintf("x*=%.4f\n",x(end))

x=linspace(-2,5,100);




