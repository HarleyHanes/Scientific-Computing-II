function [Data] = ExtractData(Data,Image)
%UNTITLED16 Summary of this function goes here
%   Detailed explanation goes here
Data.y=NaN(length(Data.x),1);
for i=1:length(Data.x)
    Data.y(i)=Image(Data.x(i,1),Data.x(i,2));
end
end

