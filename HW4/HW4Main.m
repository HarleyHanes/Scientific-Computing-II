clear; syms kappa
K=logspace(0,5);
Conv=(sqrt(K)-1)./(sqrt(K)+1);
semilogx(K,Conv)
xlabel('\kappa')
ylabel('$\frac{\sqrt{\kappa}-1}{\sqrt{\kappa}+1}$','Interpreter','latex','fontsize',13)
set(get(gca,'ylabel'),'rotation',0)
tic;

[xChol,iterChol,timeChol]=HW4PCG(A,b,tol,x0,itermax,'Cholesky');
[xSOR,iterSOR,timeSOR]=HW4PCG(A,b,tol,x0,itermax,'Cholesky');
time=toc;
