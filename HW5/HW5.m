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
%% Problem 1
clc; clear;
itermax=100; tol=1e-8;
    %Part A: General Random Matrix
        %Power Method algorithm
        A=GetAMatrix('general',15);
        q0=rand(size(A,1),1); q0=q0/norm(q0);
        [lambdaGen,~] = PowerMethod(A,q0,itermax,tol);
        %Check True Solution
        %LambdaTrue=max(eigA);
%         if abs(LambdaTrue-lambdaGen(end))>=tol*10
%             warning("%.2e difference between lambdaTrue and lambda%i",abs(LambdaTrue-lambdaGen(end)),length(lambda1))
%         end
        %Plot output
        subplot(1,4,1)
        GenErr=log(abs(lambdaGen-lambdaGen(end)));
        plot(0:length(lambdaGen)-1,GenErr)
        ylabel("log_{10}(|\lambda_{end}-\lambda_{k}|)")
        xlabel("k")
        title("Fully Random Matrix")
    %Part B: diag-dominant and symmetric
        %Diag-dominant
        A=GetAMatrix('DiagDom',15);
        q0=rand(size(A,1),1); q0=q0/norm(q0);
        [lambdaDD,~] = PowerMethod(A,q0,itermax,tol);
        %Check True Solution
        %LambdaTrue=max(eigA);
%         if abs(LambdaTrue-lambdaDD(end))>=tol*10
%             warning("%.2e difference between lambdaTrue and lambda%i",abs(LambdaTrue-lambdaGen(end)),length(lambda1))
%         end
        DDErr=log(abs(lambdaDD-lambdaDD(end)));
        subplot(1,4,2)
        plot(0:length(lambdaDD)-1,DDErr)
        xlabel("k")
        title("Diagonally Dominant Matrix")
        
        %Symmetric
        A=GetAMatrix('symmetric',15);
        q0=rand(size(A,1),1); q0=q0/norm(q0);
        [lambdaSym,~] = PowerMethod(A,q0,itermax,tol);
        %Check True Solution
        %LambdaTrue=max(eigA);
%         if abs(LambdaTrue-lambdaSym(end))>=tol*10
%             warning("%.2e difference between lambdaTrue and lambda%i",abs(LambdaTrue-lambdaGen(end)),length(lambda1))
%         end
        SymErr=log(abs(lambdaSym-lambdaSym(end)));
        subplot(1,4,3)
        plot(0:length(lambdaSym)-1,SymErr)
        xlabel("k")
        title("Symmetric Matrix")
        
        %Complex-Conjugate
        A=GetAMatrix('ImagEig',15);
        q0=rand(size(A,1),1); q0=q0/norm(q0);
        [lambdaCom,x1] = PowerMethod(A,q0,itermax,tol);
        %Check True Solution
        %LambdaTrue=max(eigA);
%         if abs(LambdaTrue-lambdaSym(end))>=tol*10
%             warning("%.2e difference between lambdaTrue and lambda%i",abs(LambdaTrue-lambdaGen(end)),length(lambda1))
%         end
        ComErr=log(abs(lambdaCom-lambdaCom(end)));
        subplot(1,4,4)
        plot(0:length(lambdaCom)-1,ComErr)
        xlabel("k")
        title("Complex Conjugates")
    
    %% Part D
        clear; itermax=1000; tol=1e-10;
        A=GetAMatrix('general',15);
        q0=rand(size(A,1),1); q0=q0/norm(q0);
        [lambda1,x1] = PowerMethod(A,q0,itermax,tol);
        [X,D,V]=eig(A);
        x1True=X(:,1);
            %Check x1 eigenvector
            fprintf("x1 is different from the true eigenvector on the order of %.2e\n",norm(x1-x1True))
        v0=rand(size(A,1),1);
        x1=x1True;
        z0=v0-(v0'*x1)*x1;        %Checked: z0 is orthoginal to x1 at tol=1e-10
        z0=z0/norm(z0);
            fprintf("z0 is orthoginal to x1 on the order of %.2e\n",z0'*x1)
        [lambda2,x2]=PowerMethod(A,z0,itermax,tol);
            %Outputting lambda2=lambda1

%% Problem 2
    %% Part a: Check Function finds lambda1
    clear; tol=1e-8; itermax=100;
    for i=1:3
        if i==1
            A=GetAMatrix('general',15);
        elseif i==2
            A=GetAMatrix('DiagDom',15);
        elseif i==3
            A=GetAMatrix('symmetric',15);
        end
        q0=rand(size(A,1),1); q0=q0/norm(q0);
        [lambdaP,~]=PowerMethod(A,q0,itermax,tol);
        [lambdaI,q,~] = InversePowerMethod(A,q0,lambdaP(end)+1e-5,itermax,tol,1);
        lambdaP=lambdaP(end);
        lambdaI=lambdaI(end);
        lambdaTrue=max(eig(A));
        if abs(lambdaTrue-lambdaI) >=tol
            warning("InversePowerMethod Unsucessful")
        end
    end
    %% Part b: Find Real eigenvalues
    clear; tol=1e-10; itermax=200;
    for j=3:-1:1
        if j==1
            A=GetAMatrix('general',15);
        elseif j==2
            A=GetAMatrix('DiagDom',15);
        elseif j==3
            A=GetAMatrix('symmetric',15);
        end
        q0=rand(size(A,1),1);q0=q0/norm(q0);
        %seperate real and imaginary eigenvalues
        RealEig=0;ImagEig=0; RealEig(1)=[];ImagEig(1)=[];
        eigA=eig(A);
        for i=1:length(A)
            if isreal(eigA(i))
                RealEig(end+1)=eigA(i);
            else
                ImagEig(end+1)=eigA(i);
            end
        end
        %Create Series of mu values along which to search
        mu=linspace(min(RealEig),max(RealEig),400)+.01;
        lambda=0; lambda(1)=[];
        for i=1:length(mu)
            [lambdaTest,~,success] = InversePowerMethod(A,q0,mu(i),itermax,tol,0);
            lambdaTest=lambdaTest(end);
            if sum(abs(lambda-lambdaTest)>=tol*10000)==length(lambda) && success  %Test if lambdaP new
                lambda(end+1)=lambdaTest;     %If new, append lambdaP
            end
        end
        lambda=sort(lambda);
        fprintf("%i real eigenvalues out of %i found\n",length(lambda),length(RealEig));
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
%clear; clc;
nsamples=50; size=20; itermax=1000;
MeanErr=zeros(3,itermax);
for i=1:nsamples
    Asym=GetAMatrix('symmetric',size);
    ASPD=GetAMatrix('SPD',size);
    ATri=GetAMatrix('SymTriDiag',size);
    [eigSym,~]=QRalgorithm(Asym,itermax);
    [eigSPD,~]=QRalgorithm(ASPD,itermax);
    [eigTri,~]=QRalgorithm(ATri,itermax);
    errSym=sum(abs((sort(eigSym)-sort(eig(Asym)))))/norm(eig(Asym));
    errSPD=sum(abs((sort(eigSPD)-sort(eig(ASPD)))))/norm(eig(ASPD));
    errTri=sum(abs((sort(eigTri)-sort(eig(ATri)))))/norm(eig(ATri));
    MeanErr(1,:)=((i-1)*MeanErr(1,:)+errSym)/i;
    MeanErr(2,:)=((i-1)*MeanErr(2,:)+errSPD)/i;
    MeanErr(3,:)=((i-1)*MeanErr(3,:)+errTri)/i;
end
semilogy(1:itermax,MeanErr)
 legend('Symmetric','SPD','Symmetric Tridiagonal')
 ylabel('||\lambda_{QR}-\lambda_{True}||_1')
 xlabel('Iterations')



