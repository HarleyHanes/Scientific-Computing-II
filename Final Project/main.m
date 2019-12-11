%Main script for final project on RBFs
clc; clear; close all;
    %BasisFunc outputs a row vector with the f(x) values at x for each xi
    
%% Simulated Data
    %Data generation
        %We will test our RBF methods on data Halton Data Sets described in RBF.pdf
        f=@(x,y)sinc(pi.*x).*sinc(pi.*y);
        dim=[0, 0; 1 1];
        numpoints=500;
        Data=GenerateHalton(f,dim,numpoints);
    
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
        
        keyboard
%%Simulation 3: Error Comparison
   
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
      
      
%         numRef=.1*numpoints;
%         %Get Reference Points
%         xiApprox=[0 0 1 1 Data.x(end-numRef:end,1)';
%                   1 0 1 0 Data.x(end-numRef:end,2)'];
%         %Plot Data and reference points
%         figure
%         plot(Data.x(:,1),Data.x(:,2),'*','MarkerSize',.2)
%         hold on
%         plot(xiApprox(1,:),xiApprox(2,:),'rs','MarkerSize',3)
%         legend('Halton Sample Points','RBF Reference Points')
%         title('Sample and Reference Points for 89 Halton Points with Corners')
% 
%         %Get Interp planar coeffecients
%         [betaPlanarInterp,~,tcost]=SolveRBFPlanarInterp(Data,basisFunc,0);
% 
%         %Get Approx coeffecients
%         betaApprox=SolveRBFApprox(Data,xiApprox,basisFunc);
% 
%         %Get Approx with Plane coeffecients 
%         betaPlanarApprox=SolveRBFPlanarApprox(Data,xiApprox,basisFunc);
% 
%         FPlanarInterp=@(x)dot([basisFunc(x,Data.x') x 1],betaPlanarInterp);
% 
%         FApprox=@(x)dot(basisFunc(x,xiApprox),betaApprox);
% 
%         FPlanarApprox=@(x)dot([basisFunc(x,xiApprox) x 1],betaPlanarApprox);
% 
% 
%         %Test on Rectinilar Grid
%         gridpoints=100;
%         x=linspace(dim(1,1),dim(2,1),gridpoints);
%         y=linspace(dim(1,2),dim(2,2),gridpoints);
%         fMesh=NaN(gridpoints);
%         InterpPlanarMesh=NaN(gridpoints);
%         ApproxMesh=NaN(gridpoints);
%         ApproxPlanarMesh=NaN(gridpoints);
%         for i=1:gridpoints
%             for j=1:gridpoints
%                 fMesh(i,j)=f(x(i),y(j));
%                 InterpPlanarMesh(i,j)=FPlanarInterp([x(i) y(j)]);
%                 ApproxMesh(i,j)=FApprox([x(i) y(j)]);
%                 ApproxPlanarMesh(i,j)=FPlanarApprox([x(i) y(j)]);
%             end
%         end
%         %Compute Error
%         ErrApprox=sum(sum(sqrt((ApproxMesh-fMesh).^2)))/(gridpoints^2);
%         ErrPlanarApprox=sum(sum(sqrt((ApproxPlanarMesh-fMesh).^2)))/(gridpoints^2);
% 
%         figure
%         mesh(x,y,InterpPlanarMesh)
% 
%         figure
%         mesh(x,y,fMesh)
%         title('$f(x)=sinc(\pi x)*sinc(\pi y)$','Interpreter','LaTex')
% 
%         figure
%         mesh(x,y,ApproxMesh)
%         titleApprox=sprintf('RBF Approx., Mean err=%.3e',ErrApprox);
%         title(titleApprox)
% 
%         figure
%         mesh(x,y,ApproxPlanarMesh)
%         titlePlanarApprox=sprintf('RBF Planar Approx., Mean err=%.3e',ErrPlanarApprox);
%         title(titlePlanarApprox)
% 
%         %Asses Sample Point effect 
% 
%         clear ErrApprox ErrPlanarApprox betaApprox betaPlanarApprox xiApprox;
%         keyboard
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

    [betaApprox,condApprox(iRef)]=SolveRBFApprox(Data,xiApprox,basisFunc);

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
alpha=.01;
basisFunc=@(x,xi)exp(-(alpha.*sqrt(sum((x'-xi).^2))).^2);
%LoadImage
ImageTrue=im2double(rgb2gray(imread('../../StreetCar.jpg')));
subplot(1,3,1)
imshow(ImageTrue)
title('True Image')

imagDim=size(ImageTrue);

%DistortImage
%1)Add Noise
ImageNoisy=ImageTrue+(rand(size(ImageTrue))-.5)*mean(mean(ImageTrue)/1.5);
subplot(1,3,2)
imshow(ImageNoisy)
title('Image with Noise')

%2)AddCreases
ImageCreased=ImageTrue;
for i=100:100:size(ImageTrue,1)-100
    ImageCreased(i-5:i+5,:)=ones(11,size(ImageTrue,2));
end
for j=100:100:size(ImageTrue,2)-100
    ImageCreased(:,i-5:i+5)=ones(size(ImageTrue,1),11);
end
subplot(1,3,3)
imshow(ImageCreased)
title('Image with Creases')

ImageData=ImageToData(ImageTrue);

%Generate Reference Points
% NumRef=100;
% DimRatio=imagDim(1)/imagDim(2);
numXxi=90;numYxi=70;
xix=round(linspace(1,imagDim(1),numXxi+2));
xiy=round(linspace(1,imagDim(2),numYxi+2));
xix=xix(2:end-1);xiy=xiy(2:end-1);
xiImage=[0;1];
xiImage(:,1)=[];
for i=1:length(xix)
    xiImage=[xiImage [xix(i).*ones(1,length(xiy));xiy]];
end

%Translate Images to Data struct
% RawDataTrue=ImageToData(ImageTrue);
% RawDataNoisy=ImageToData(ImageNoisy);
% RawDataCreased=ImageToData(ImageCreased);

%Create Halton points
pointstep=5;
PointIndices=1:pointstep:imagDim(1)*imagDim(2);
for i=1:length(PointIndices)
    DataTrue.x(i,:)=ImageData.x(PointIndices(i),:);
end
% p=haltonset(2);
% DataTrue.x=ceil(net(p,numpoints).*imagDim);
% DataTrue.x(1,:)=[];
DataNoisy.x=DataTrue.x;
DataCreased.x=DataTrue.x;

DataTrue=ExtractData(DataTrue,ImageTrue);
DataNoisy=ExtractData(DataNoisy,ImageNoisy);
DataCreased=ExtractData(DataCreased,ImageCreased);

%True Fitting
    %Get Approx coeffecients
    %betaTrueApprox=SolveRBFApprox(DataTrue,xiImage,basisFunc);
    %Get Approx with Plane coeffecients 
    betaTruePlanarApprox=SolveRBFPlanarApprox(DataTrue,xiImage,basisFunc);

%Noisy Fitting
    %Get Approx coeffecients
    %betaNoisyApprox=SolveRBFApprox(DataNoisy,xiImage,basisFunc);
    %Get Approx with Plane coeffecients 
    %betaNoisyPlanarApprox=SolveRBFPlanarApprox(DataNoisy,xiImage,basisFunc);

%Creased Fitting
    %Get Approx coeffecients
    %betaCreasedApprox=SolveRBFApprox(DataCreased,xiImage,basisFunc);
    %Get Approx with Plane coeffecients 
   % betaCreasedPlanarApprox=SolveRBFPlanarApprox(DataCreased,xiImage,basisFunc);



FTrueApprox=@(x)dot(basisFunc(x,xiImage),betaTrueApprox);
FTruePlanarApprox=@(x)dot([basisFunc(x,xiImage) x 1],betaTruePlanarApprox);

FNoisyApprox=@(x)dot(basisFunc(x,xiImage),betaNoisyApprox);
FNoisyPlanarApprox=@(x)dot([basisFunc(x,xiImage) x 1],betaNoisyPlanarApprox);

FCreasedApprox=@(x)dot(basisFunc(x,xiImage),betaCreasedApprox);
FCreasedPlanarApprox=@(x)dot([basisFunc(x,xiImage) x 1],betaCreasedPlanarApprox);

DataTrueApprox=ImageData;
DataTruePlanarApprox=ImageData;

for i=1:length(ImageData.x)
    %DataTrueApprox.y(i)=FTrueApprox(ImageData.x(i,:));
    DataTruePlanarApprox.y(i)=FTruePlanarApprox(ImageData.x(i,:));
end
    
%TrueApprox=DataToImage(DataTrueApprox);
TruePlanarApprox=DataToImage(DataTruePlanarApprox);
% figure
% imshow(TrueApprox)
figure
imshow(TruePlanarApprox)
% 
% DataCreasedApprox=FCreasedApprox(ImageData.x);
% DataCreasedPlanarApprox=FCreasedApprox(ImageData.x);
% 
% DataNoisyApprox=FNoisyApprox(ImageData.x);
% DataNoisyPlanarApprox=FNoisyApprox(ImageData.x);

