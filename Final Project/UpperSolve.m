function xU = UpperSolve(U,b)
sizecheck=size(U);
if sizecheck(1)~=sizecheck(2)
    fprintf("Error!! Entered matrix is not square")
    keyboard
end
n=sizecheck(1);
xU=NaN(size(b));
for i=n:-1:1
    if i==n
        xU(i)=b(i)/U(i,i);
    else
        j=n:-1:i+1;
        xU(i)=(b(i)-U(i,j)*xU(j))/U(i,i);
    end
end
end