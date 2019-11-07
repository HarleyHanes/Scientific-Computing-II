%Homework 5- Scientific Computing II
%Harley Hanes
clear; close all;clc;
%% Matrix Creation and Settings
%Construct A
 A=GetAMatrix('general',15);
[X,D,Xinv]=eig(A);
eigA=sort(diag(D));
RealEig=0;ImagEig=0;
for i=1:length(eigA)
    if isreal(eigA(i))
        RealEig(end+1)=eigA(i);
    else
        ImagEig(end+1)=eigA(i);
    end
end
RealEig(1)=[];ImagEig(1)=[];
LambdaTrue=max(eigA);
%Set iteration settings
itermax=100;
tol=1e-10;
%% Problem 1
    %Compute Lambda1
        %Power Method algorithm
        q0=rand(size(A,1),1); q0=q0/norm(q0);
        [lambda1,x1] = PowerMethod(A,q0,itermax,tol);
        %Check True Solution
        if abs(LambdaTrue-lambda1(end))>=tol*10
            warning("%.2e difference between lambdaTrue and lambda%i",abs(LambdaTrue-lambda1(end)),length(lambda1))
        end

        %Plot output
        LambdaErr=log(abs(lambda1-lambda1(end)));
        plot(0:length(lambda1)-1,LambdaErr)
        ylabel("log_{10}(|\lambda_{end}-\lambda_{k}|)")
        xlabel("k")
        title("Power Method")
    
    %ComputeLambda2-Currently just outputting lambda1
        %x1=Xinv(1,:)';
        %x1=X(:,1);
        v0=rand(size(A,1),1);
        z0=v0-(v0'*x1)*x1;        %Checked: z0 is orthoginal to x1 at tol=1e-10
        z0=z0/norm(z0);
        [lambda2,x2]=PowerMethod(A,z0,itermax,tol);
            %Outputting lambda2=lambda1

%% Problem 2
    %% Part a: Check Function finds lambda1 for mu=0
    q0=rand(size(A,1),1);
    [lambdaP,q,~] = InversePowerMethod(A,q0,lambda1(end)+1e-5,itermax,tol,1);
    lambdaP=lambdaP(end);
    if abs(lambdaP-lambda1(end)) >= tol
        warning("InversePowerMethod Unsucessful for u=0")
    end
    %% Part b: Find Real eigenvalues
    mu=linspace(lambda1(end),-lambda1(end),400);
    lambdaReal=lambda1(end);
    for i=1:length(mu)
        [lambdaTest,~,success] = InversePowerMethod(A,q0,mu(i),itermax,tol,0);
        if sum(abs(lambdaReal-lambdaTest(end))>=tol*100)==length(lambdaReal) && success  %Test if lambdaP new
            lambdaReal(end+1)=lambdaTest(end);     %If new, append lambdaP
        end
    end
    lambdaReal=sort(lambdaReal);
    fprintf("%i real eigenvalues out of %i found\n",length(lambdaReal),length(RealEig));
    %This needs better implementation to account for when not all real
    %eigenvalues are found
    if sum(lambdaReal-RealEig>=tol*100)
        warning("Real eigenvalues found not matching true real eigenvalues")
        disp(lambdaReal-RealEig)
    end
    
    
    %% Part c: Find lambda1 and lambda 2 given that they're complex
    %conjugates
    %Get A matrix
    A=GetAMatrix('ImagEig',15);
    q0=rand(15,1); q0=q0/norm(q0);
    %Use PowerMethod to get guess of the real part
    [RealGuess,~] = PowerMethod(A,q0,10,tol);
    %Get mu closer to positive complex conjugate
     mu=RealGuess(end)+2i;
     [lambda1,~]=InversePowerMethod(A,q0,mu,itermax,tol,1);
     lambda1=lambda1(end);       %Ignore previous entries
     %Select mu very close to lambda2
     mu=lambda1-(2*imag(lambda1)-.1)*1i;
     [lambda2,~]=InversePowerMethod(A,q0,mu,itermax,tol,1);
     lambdaTrue=sort(eig(A));
     fprintf("Error in lambda1: %.3e\n",lambdaTrue(end)-lambda1(end))
     fprintf("Error in lambda2: %.3e\n",lambdaTrue(end-1)-lambda2(end))
     
%% Problem 3
clear;
SizeVec=10:10:100;
tol=1e-6; itermax=10000;
Eig=NaN(3,length(SizeVec));
for i=1:length(SizeVec)
    size=SizeVec(i);
    Asym=GetAMatrix('symmetric',size);
    ASPD=GetAMatrix('SPD',size);
    ATri=GetAMatrix('SymTriDiag',size);
    [eigSym,iterSym]=QRalgorithm(Asym,itermax,tol);
    [eigSPD,iterSPD]=QRalgorithm(ASPD,itermax,tol);
    [eigTri,iterTri]=QRalgorithm(ATri,itermax,tol);
    Eig(:,i)=[iterSym iterSPD iterTri]';
end
plot(SizeVec,Eig)
legend('Symmetric','SPD','Symmetric Tridiagonal')

%Results
%Only works on symmetric matrices, faster on symmetric than SPD

