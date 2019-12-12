function [Data] = ExtractData(Data,Image)
%ExtractData Given a Data.x datapoints and Image matrix, extracts Data.y
Data.y=NaN(length(Data.x),1);
for i=1:length(Data.x)
    Data.y(i)=Image(Data.x(i,1),Data.x(i,2));
end
end

