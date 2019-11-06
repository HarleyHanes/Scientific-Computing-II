function Anew=GetImagEigA(size) %Find an A with a real maximum eigenvalue
%Find random A
A=rand(size);
[X,D,~]=eig(A);
%find max eigenvalue
maxEig=max(abs(diag(D)));
if ~isreal(maxEig)
    Anew=A;
else
    realComponent=(rand(1)+2)*maxEig;  %We know the real component larger than maxEig
    imagComponent=(rand(1)+2);
    lambda1=realComponent+imagComponent*1i;   %|lambda1|>maxEig
    lambda2=realComponent-imagComponent*1i;   %|lambda2|>maxEig
    DD=D;
    for i=1:size
        if ~isreal(D(i,i))  %Find first imaginary eigenvalue
            DD(i,i)=lambda1;
            DD(i+1,i+1)=lambda2;
            break           %Stop searching after one complex eigenvalue
        end
    end
    Anew=X*DD/X;
    Anew=real(Anew);      %Extract real parts to eliminate roudnoff errors
end
%Check Anew has correct eigenvalues
Dnew=sort(eig(Anew));
if isreal(Dnew(end)) || ~isreal(Dnew(end)*Dnew(end-1))
    error("Anew does not have imaginary maximum eigenvalues")
end
end