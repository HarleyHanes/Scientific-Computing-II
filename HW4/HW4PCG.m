function [x,iter,err,time,Conv] = HW4PCG(A,b,x0,tol,itermax,Preconditioner)
%HW4PCG Preconditioned Conjugate Gradient for HW4 of Scientifc Computing II
%Authors: Harley Hanes, Elliot Hill-10/21/19
%References: Code based on examples from Barlow J. 2019. CSE/Math 456: Lecture 
    %20. The Pennsylvania State University. http://www.cse.psu.edu/~b58/cse456/lecture20.pdf
%   Outputs
%       x: Final Solution of Ax=b found under PCG
%       iter: number of iterations to reach x
%   Inputs
%       A: nxn SPD matrix
%       b: nx1 vector
%       x0: nx1 vector, initial guess for x
%       tol: scalar for tolerance of function
%       itermax: scalar for maximum number of iterations run
%       Preconditioner: string selecting preconditioning method
%           Supported Preconditioners: Incomplete Cholesky

%% Check Inputs
%A Checks
[m,n]=size(A);
    %Check A Square
    if m~=n
        warning("Error!! A non-square")
        keyboard
    end
    %Check A SPD
    if (A')~=A
        warning("A is not symmetric")
        keyboard
    end
    EigA=eig(A'*A);
    NegEigA=EigA(EigA<0);
    if isempty(NegEigA)==0
        NegEigAList=sprintf("%.3e  ",NegEigA);
        warning("Error!! A has negative eigenvalue(s): %s",NegEigAList)
        keyboard
    end
%b Checks
blength=length(b);
    if blength~=m
        warning("A and b sizes don't match")
        keyboard
    end
%x0 Checks
x0length=length(x0);
    if x0length~=m
        warning("A and x0 sizes don't match")
        keyboard
    end
%% Apply Preconditioner
switch Preconditioner
    case "Cholesky"
        %Check A Irreducibly Diag Dominant
            DiagDomFlag=0;
            StrictDomFlag=0;
            for i=1:length(A)
                if DiagDomFlag==0 && 2*abs(A(i,i)) < sum(abs(A(:,i)))
                    DiagDomFlag=1;
                elseif StrictDomFlag==1 && 2*abs(A(i,i)) > sum(abs(A(:,i)))
                    StrictDomFlag=0;
                end
            end
            if StrictDomFlag==1 || DiagDomFlag==1
                warning('A is not irreducibly diagonally dominant')
                keyboard
            end
        Asparse=sparse(A);
        %L=iCholDecomp(A);
        L=ichol(Asparse,struct('droptol',10e-15,'type','nofill'));
        %Linv=inv(L);
        M=L'*L;
        %Minv=Linv*Linv';
    case "Jacobi"
        L=diag(sqrt(diag(A)));
        M=diag(diag(A));
    case "None"
        M=1;
        L=1;
    otherwise 
        warning("Error!! Preconditioner type not recognized")
end
%kappa=cond(L\A/L');
kappa=cond(inv(L)*A*inv(L'));
Conv=(sqrt(kappa)-1)/(sqrt(kappa)+1);

%% Conjugate Gradient
%[x,err,iter,~]=conjgrad(A,b,x0,M,itermax,tol);
%intialize variables
    iter=1;
     r=b-A*x0; x=x0;
tic;
 while 1  
   z=M\r; rho=r'*z;
    if iter>1
        beta=rho/rho1;
        p=z+beta*p;
    else
        p=z;
    end
    q=A*p;
    alpha=rho/(p'*q);
    x=x+alpha*p;
    r=r-alpha*q;
    rho1=rho;
    err=norm(r);
    if err <= tol %Successful convergence
        break
    end
    if iter == itermax   %UnSuccessful convergence
        warning("PCG with %s preconditioner did not converge to %.3e in %i iterations",Preconditioner,tol,itermax)
        break
    end
    if mod(iter,10)==0
        fprintf('Iteration: %i\n',iter);
        fprintf('Norm of Residual: %e\n',err);
    end 
    iter=iter+1;
 end
time=toc;
end
