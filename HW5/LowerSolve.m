function xL = LowerSolve(L,b)
sizecheck=size(L);
if sizecheck(1)~=sizecheck(2)
    fprintf("Error!! Entered matrix is not square")
    keyboard
end
n=sizecheck(1);
xL=NaN(size(b));
for i=1:n
    if i==1
        xL(i)=b(i)/L(i,i);
    else
        j=1:i-1;
        xL(i)=(b(i)-L(i,j)*xL(j))/L(i,i);
    end
end
end