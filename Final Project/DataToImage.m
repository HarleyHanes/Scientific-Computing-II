function [Image] = DataToImage(Data)
%DataToImage Takes a Data struct and Turns it to image matrix
%   not dependent on ordering of Data points but won't check Data is full
%Image=NaN(dim);
for i=1:length(Data.y)
%     x1index=mod(i,dim(1));
%     x2index=floor(i/dim(1));
    Image(Data.x(i,1),Data.x(i,2))=Data.y(i);
    Data.y(i)=NaN;
end
if sum(isnan(Data.y))~=length(Data.y)
    warning('%i indexes not taken of Data.y',length(Data.y)-sum(isnan(Data.y)))
end

