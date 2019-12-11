function [Data] = ImageToData(Image)
%ImagToScatter Takes an Image matrix and converts it to Dat struct
%   Sorts Image by columns then rows (ie finishes all row 1 then does row 2)
Data.x=NaN(size(Image,1)*size(Image,2),2);
Data.y=NaN(size(Image,1)*size(Image,2),1);
for i=1:size(Image,1)
    Data.x((i-1)*size(Image,2)+1:i*size(Image,2),1)=i*ones(size(Image,2),1);
end

for i=1:size(Image,1)
    Data.x((i-1)*size(Image,2)+1:i*size(Image,2),2)=1:size(Image,2);
end
   
for ix2=1:size(Image,2)
    for ix1=1:size(Image,1)
        Data.y(ix2+size(Image,2)*(ix1-1))=Image(ix1,ix2);
    end
end
%CheckOuputs
if length(Data.x)~=size(Image,1)*size(Image,2)
    warning('Image has %i points, x has %i',size(Image,1)*size(Image,2),length(Data.x))
end
if length(Data.y)~=size(Image,1)*size(Image,2)
    warning('Image has %i points, y has %i',size(Image,1)*size(Image,2),length(Data.y))
end
end

