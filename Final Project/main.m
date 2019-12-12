%Main script for final project on RBFs
clc; clear; close all;
%% Simulated Data
    %Data generation
        %We will test our RBF methods on data Halton Data Sets described in RBF.pdf
        f=@(x,y)sinc(pi.*x).*sinc(pi.*y);
        dim=[0, 0; 1 1];
        numpoints=500;
        Data=GenerateHalton(f,dim,numpoints);
        plot(Data.x(:,1),Data.x(:,2),'*','MarkerFaceColor','b')
        hold on
        plot(Data.x(48:72,1),Data.x(48:72,2),'o','MarkerSize',10,'MarkerFaceColor','r')
    
    %% Simulation 1- Stability of Interpolations
        %Test Stability of no Regularization vs Regularization at different
        %alpha
         gridpoints=[100 100];
%         x=linspace(dim(1,1),dim(2,1),gridpoints);
%         y=linspace(dim(1,2),dim(2,2),gridpoints);
%         fMesh=NaN(gridpoints);
%         for i=1:gridpoints
%             for j=1:gridpoints
%                 fMesh(i,j)=f(x(i),y(j));
%             end
%         end
         alphaInterp=.25:.25:15;
        maxErr=NaN(length(alphaInterp),2);
        for ialpha=1:length(alphaInterp)
            interpMesh=NaN(gridpoints);
            interpRegMesh=NaN(gridpoints);
            alpha=alphaInterp(ialpha);
            basisFunc=@(x,xi)exp(-(alpha.*sqrt(sum((x'-xi).^2))).^2);
            [beta,cond(ialpha,1),tcost(ialpha,1),success(ialpha,1)]=SolveRBFPlanarInterp(Data,basisFunc,0);
            [betaReg,cond(ialpha,2),tcost(ialpha,2),success(ialpha,2)]=SolveRBFPlanarInterp(Data,basisFunc,1);
            if success(ialpha,1)==1
                FInterp=@(x)dot([basisFunc(x,Data.x') x 1],beta); 
                maxErr(ialpha,1) = EstimateError(gridpoints,dim,f,FInterp,'Max');
            end
            if success(ialpha,2)==1
                FInterpReg=@(x)dot([basisFunc(x,Data.x') x 1],betaReg);
                maxErr(ialpha,2) = EstimateError(gridpoints,dim,f,FInterpReg,'Max');
            end
            
            
        end
        turnpoint=sum(success(:,1)==0);
        figure
        semilogy(alphaInterp,cond)
        hold on
        %semilogy(alphaInterp(turnpoint),cond(turnpoint,1),'*')
        xlabel('\alpha')
        ylabel('Condition Number')
        
        
        figure
        plot(alphaInterp,tcost(:,2))
        hold on
        plot(alphaInterp(turnpoint:end),tcost(turnpoint:end,1))
        xlabel('\alpha')
        ylabel('Computation Time (s)')
        legend('Un-Regularized','Regularized')
        
        figure
        semilogy(alphaInterp,maxErr(:,2))
        hold on
        semilogy(alphaInterp(turnpoint:end),maxErr(turnpoint:end,1))
        xlabel('\alpha')
        ylabel('Max Error')
  
        
        
    %% Simulation 2-Comparison of computation time and MSE
    %Gauss basis with L2 norm, x is 1x2 while xi is 2xM
    clear; clc;
        
        gridpoints=[100 100];
        f=@(x,y)sinc(pi.*x).*sinc(pi.*y);
        dim=[0, 0; 1 1];
        numpoints=500;
        Data=GenerateHalton(f,dim,numpoints);
        alphaApprox=1;
        basisFunc=@(x,xi)exp(-(alphaApprox.*sqrt(sum((x'-xi).^2))).^2);
        numRef=5:20:450;
        cond=NaN(length(numRef),1); tcost=NaN(length(numRef),1);
        for iRef=1:length(numRef)
            Ref=numRef(iRef);
            xiApprox=[0 0 1 1 Data.x(end-Ref:end,1)';
                      1 0 1 0 Data.x(end-Ref:end,2)'];
            [beta,cond(iRef),tcost(iRef),success(iRef)]=SolveRBFPlanarApprox(Data,xiApprox,basisFunc);
            if success(iRef)==1
                FApprox=@(x)dot([basisFunc(x,xiApprox) x 1],beta); 
                maxErr(iRef) = EstimateError(gridpoints,dim,f,FApprox,'Max');
            end
        end
        figure
        semilogy(numRef,cond)
        xlabel('Number of Reference Points')
        ylabel('Condition Number')
        
        
        figure
        plot(numRef,tcost)
        xlabel('Number of Reference Points')
        ylabel('Computation Time (s)')
        
        figure
        plot(numRef,maxErr)
        xlabel('Number of Reference Points')
        ylabel('Max Error')
     
%% Simulation 3: Error Comparison
        clear; clc;
        gridpoints=[100 100];
        f=@(x,y)sinc(pi.*x).*sinc(pi.*y);
        dim=[0, 0; 1 1];
        numpoints=500;
        Data=GenerateHalton(f,dim,numpoints);
  
        
        alphaApprox=1;
        alphaInterp=3;
        numRef=50;
        
        basisFuncApprox=@(x,xi)exp(-(alphaApprox.*sqrt(sum((x'-xi).^2))).^2);
        basisFuncInterp=@(x,xi)exp(-(alphaInterp.*sqrt(sum((x'-xi).^2))).^2);

        xiApprox=[0 0 1 1 Data.x(end-numRef:end,1)';
          1 0 1 0 Data.x(end-numRef:end,2)'];
         %Get Interp Err
         [betaPlanarInterp,~,tcost]=SolveRBFPlanarInterp(Data,basisFuncInterp,1);
         FPlanarInterp=@(x)dot([basisFuncInterp(x,Data.x') x 1],betaPlanarInterp);
         ErrInterp=EstimateError(gridpoints,dim,f,FPlanarInterp,'all');
      
         %Get Approx with Plane coeffecients 
         betaPlanarApprox=SolveRBFPlanarApprox(Data,xiApprox,basisFuncApprox);
         FPlanarApprox=@(x)dot([basisFuncApprox(x,xiApprox) x 1],betaPlanarApprox);
         ErrApprox=EstimateError(gridpoints,dim,f,FPlanarApprox,'all');
         
         %Make Histograms
         figure
         histogram(ErrInterp/length(ErrInterp))
         xlabel('Interpolation Estimate Error')
         
         figure
         histogram(ErrApprox/length(ErrApprox))
         xlabel('Approximation Estimate Error')
         
         %Plot Curves
        gridpoints=100;
        x=linspace(dim(1,1),dim(2,1),gridpoints);
        y=linspace(dim(1,2),dim(2,2),gridpoints);
        fMesh=NaN(gridpoints);
        InterpPlanarMesh=NaN(gridpoints);
        ApproxMesh=NaN(gridpoints);
        ApproxPlanarMesh=NaN(gridpoints);
        for i=1:gridpoints
            for j=1:gridpoints
                fMesh(i,j)=f(x(i),y(j));
                InterpPlanarMesh(i,j)=FPlanarInterp([x(i) y(j)]);
                ApproxPlanarMesh(i,j)=FPlanarApprox([x(i) y(j)]);
            end
        end
figure
mesh(x,y,fMesh)
title('$f(x)=sinc(\pi x)*since(\pi y)$','Interpreter','Latex')

figure
mesh(x,y,InterpPlanarMesh)
title('RBF Interpolation')

figure
mesh(x,y,ApproxPlanarMesh)
title('RBF Approximation')
        
        
        
%%

%Set Number of reference points
numRef=(5:10:length(Data.y))';

%Allocate Error Matrices
ErrApprox=NaN(length(numRef),3);
ErrPlanarApprox=NaN(length(numRef),3);

%Allocate Condition number vectors
condApprox=NaN(length(numRef),1);
condPlanarApprox=NaN(length(numRef),1);

%Allocate t vectors
tApprox=NaN(length(numRef),1);
tPlanarApprox=NaN(length(numRef),2);


for iRef=1:length(numRef)
    xiApprox=[0 0 1 1 Data.x(end-numRef(iRef)+4:end,1)';
          1 0 1 0 Data.x(end-numRef(iRef)+4:end,2)'];
      
      
    %Get Approx coeffecients

    [betaApprox,condApprox(iRef)]=SolveRBFApprox(Data,xiApprox,basisFuncApprox);

    %Get Approx with Plane coeffecients 
  
    [betaPlanarApprox,condPlanarApprox(iRef),tPlanarApprox(iRef,:)]=SolveRBFPlanarApprox(Data,xiApprox,basisFunc);

    
    %Define RBF
    FApprox=@(x)dot(basisFunc(x,xiApprox),betaApprox);
    FPlanarApprox=@(x)dot([basisFunc(x,xiApprox) x 1],betaPlanarApprox);
    
    %Solve on Grid
    for i=1:gridpoints
        for j=1:gridpoints
            ApproxMesh(i,j)=FApprox([x(i) y(j)]);
            ApproxPlanarMesh(i,j)=FPlanarApprox([x(i) y(j)]);
        end
    end
    
    %Check Error
    ErrApproxMat=sqrt((ApproxMesh-fMesh).^2);
    ErrPlanarApproxMat=sqrt((ApproxPlanarMesh-fMesh).^2);
    
    ErrApprox(iRef,1)=min(min(ErrApproxMat));
    ErrPlanarApprox(iRef,1)=min(min(ErrPlanarApproxMat));
    
    ErrApprox(iRef,2)=sum(sum(ErrApproxMat))/(gridpoints^2);
    ErrPlanarApprox(iRef,2)=sum(sum(ErrPlanarApproxMat))/(gridpoints^2);
    
    ErrApprox(iRef,3)=max(max(ErrApproxMat));
    ErrPlanarApprox(iRef,3)=max(max(ErrPlanarApproxMat));
end

%Plot Results
figure
%Error Estimates
subplot(3,2,1)
semilogy(numRef,ErrApprox)
ylabel('Err on Grid')
title('RBF Approx.')

subplot(3,2,2)
semilogy(numRef,ErrPlanarApprox)
legend('Max err','Mean err','Min err')
title('RBF Planar Approx.')

%Condition Number
subplot(3,2,3)
semilogy(numRef,condApprox)
ylabel('Cond(A^TA)')

subplot(3,2,4)
semilogy(numRef,condPlanarApprox)

%Tcost
subplot(3,2,5)
%plot(numRef,tApprox)
xlabel('# Reference Points')
ylabel('Computation time')

subplot(3,2,6)
plot(numRef,tPlanarApprox)
xlabel('# Reference Points')




%% Test 2: Image reconstructon
clc;clear;
alpha=5;
basisFunc=@(x,xi)exp(-(alpha.*sqrt(sum((x'-xi).^2))).^2);
%LoadImage
ImageTrue=im2double(rgb2gray(imread('../../StreetCar.jpg')));
figure
imshow(ImageTrue)

imagDim=size(ImageTrue);

%DistortImage
%1)Add Noise
ImageNoisy=ImageTrue+(rand(size(ImageTrue))-.5)*mean(mean(ImageTrue)*5);
figure
imshow(ImageNoisy)


ImageData=ImageToData(ImageTrue);

numXxi=70;numYxi=45;
xix=round(linspace(1,imagDim(1),numXxi+2));
xiy=round(linspace(1,imagDim(2),numYxi+2));
xix=xix(2:end-1);xiy=xiy(2:end-1);
xiImage=[0;1];
xiImage(:,1)=[];
for i=1:length(xix)
    xiImage=[xiImage [xix(i).*ones(1,length(xiy));xiy]];
end

%Translate Images to Data struct

%Create Halton points
pointstep=1;
PointIndices=1:pointstep:imagDim(1)*imagDim(2);
for i=1:length(PointIndices)
    DataTrue.x(i,:)=ImageData.x(PointIndices(i),:);
end
DataNoisy.x=DataTrue.x;
DataNoisy=ExtractData(DataNoisy,ImageNoisy);
%DataCreased=ExtractData(DataCreased,ImageCreased);

    %betaNoisyPlanarApprox=SolveRBFPlanarApprox(DataNoisy,xiImage,basisFunc);
    betaNoisyPlanarInterp=SolveRBFPlanarInterp(DataNoisy,basisFunc,1);

FNoisyPlanarInterp=@(x)dot([basisFuncInterp(x,Data.x') x 1],betaNoisyPlanarApprox);


for i=1:length(ImageData.x)
    DataNoisyPlanarApprox.y(i)=FNoisyPlanarInterp(ImageData.x(i,:));
end
DataNoisyPlanarApprox.x=ImageData.x;
    
NoisyPlanarApprox=DataToImage(DataNoisyPlanarApprox);

figure
imshow(NoisyPlanarApprox);

