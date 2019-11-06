%% Homework 3
clc; clear; close all;
%Problem Setup
%nvec=linspace(50,80,7);
tolvec=[10^(-2)]; %10^(-3) 10^(-4) 10^(-6) 10^(-7) 10^(-8)];
nvec=50;
itermax=1000;
w=1;
SolSave=cell(length(tolvec),length(nvec));
for ii=1:length(tolvec)
    for jj=1:length(nvec)
        n=nvec(jj);
        tol=tolvec(ii);
        tstart=clock;
        m=n-1;
        x0=zeros(m^2,1);
        x=linspace(0,1,n+1)';
        y=linspace(0,1,n+1);
        Sol=NaN(m,m);
        b=zeros(m^2,1);
        U=zeros(m^2);
        
        %Set boundaries
        x0bound = @(x)  0*x;
        x1bound = @(x)  .1*sin(pi*y).^2;
        y0bound=  @(x)  2.*x.*(1-x);
        y1bound = @(x)  -x.*(1-x);
        
        x0vec=x0bound(y)';
        x1vec=x1bound(y)';
        y0vec=y0bound(x);
        y1vec=y1bound(x);
%         plot(x,x0vec)
%         hold on
%         plot(x,x1vec)
%         plot(y,y0vec)
%         plot(y,y1vec)
%         legend
%         figure
        
        % Set up b
        b(1:m)=b(1:m)-y0vec(2:end-1);
        b(m*(m-1)+1:m^2)=b(m*(m-1)+1:m^2)-y1vec(2:end-1);
        for i=1:m
            b(i*m)=b(i*m)-x1vec(i+1);
            b(1+m*(i-1))=b(1+m*(i-1))-x0vec(i+1);
        end
        %b=zeros(m^2,1);
        %Set up U
        for i=1:m^2 %Position Assessing
            for j=1:m^2 %Position Dependencies
                if (mod(j,m)==0 && mod(i,m)==1)|| (mod(j,m)==1 && mod(i,m)==0)
                     U(i,j)=0; 
                elseif j==i-m || j==i-1 || j==i+1 || j==i+m
                    U(i,j)=1;
                elseif i==j
                    U(i,j)=-4;
                end

            end

        end
        %[xSol,iter]=SORAdjusted(U,b,x0,w,tol,itermax);
        [xSol,iter,err]=HW4PCG(U,b,x0,tol,itermax,'Jacobi');
        %xSol=U\b;
        for i=1:m
            Sol(i,:)=xSol((i-1)*m+1:i*m);
        end
        Sol=[zeros(1,m);Sol; zeros(1,m)];
        Sol=[zeros(n+1,1) Sol zeros(n+1,1)];
        Sol(:,1)=x0vec;
        Sol(:,end)=x1vec;
        Sol(end,2:m+1)=Sol(end,2:m+1)+y1vec(2:m+1)';
        Sol(1,2:m+1)=Sol(1,2:m+1)+y0vec(2:m+1)';
        tend=clock;
        t(ii,jj)=sum((tend(end-1:end)-tstart(end-1:end)).*[60 1]);
        itersave(ii,jj)=iter;
        if iter<itermax-1
            Success(ii,jj)=1;
        else 
            Success(ii,jj)=0;
        end
         SolSave{ii,jj}=Sol;
          surf(x,y,Sol)
          xlabel('x')
          ylabel('y')
    end
end
subplot(2,1,1)
semilogx(tolvec,t)
ylabel('Computation time')
subplot(2,1,2)
semilogx(tolvec,itersave)
ylabel('Iterations')
xlabel('Tolerance')

function [x,iter]=SOR(A,b,x0,w,tol,nmax)
[n,~]=size(A);
iter=0;xold=x0;x=NaN(length(x0),1);
r=b-A*x0;rnorm=norm(r);
    while rnorm > tol && iter <nmax
        iter=iter+1;
        for i=1:n
            s=0;
            for j=1:i-1
                s=s+A(i,j)*x(j);
            end
            for j=i+1:n
                s=s+A(i,j)*xold(j);
            end
            x(i,1)=w*(b(i)-s)/A(i,i)+(1-w)*xold(i);
        end
        xold=x;r=b-A*x; rnorm=norm(r);
    end
end
function [x,iter]=SORAdjusted(A,b,x0,w,tol,nmax)
%% Adjust SOR algorithm for structure of laplace problem to increase 
[n,~]=size(A);
iter=0;xold=x0;x=NaN(length(x0),1);
r=b-A*x0;rnorm=norm(r);
    while rnorm > tol && iter <nmax
        iter=iter+1;
        for i=1:n
            s=0;
            jL=[i-sqrt(n) i-1];
            jM=[i+1 i+sqrt(n)];
            jL=jL(jL(jL<i)>0);
            jM=jM(jM(jM<=n)>i);
            for j=jL
                %disp(j)
                s=s+A(i,j)*x(j);
            end
            for j=jM
                %disp(j)
                s=s+A(i,j)*xold(j);
            end
            x(i,1)=w*(b(i)-s)/A(i,i)+(1-w)*xold(i);
            %keyboard
        end
        xold=x;r=b-A*x; rnorm=norm(r);
    end
%computation time
%     [n,~]=size(A);
%     iter=0;xold=x0;x=zeros(length(x0),1);
%     r=b-A*x0;rnorm=norm(r);
%     while rnorm > tol && iter <nmax
%        iter=iter+1;
%         for i=1:n
%             s=0;
%             jj=[i-sqrt(n) i-1 i+1 i+sqrt(n)];
%             
%             jl=jj(jj(jj<i)>0);
%             jm=jj(jj(jj<n)>i);
%             
%             s=s+sum(x(jl));
%             s=s+sum(xold(jm));
%             x(i,1)=w*(b(i)-s)/A(i,i)+(1-w)*xold(i);
%         end
%     xold=x;r=b-A*x; rnorm=norm(r);
%     end

end
