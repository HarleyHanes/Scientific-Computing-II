function A=GetAMatrix(type,size) %Find an A with a real maximum eigenvalue
switch type
    case 'ImagEig'
        %Find random A
        A=rand(size);
        [X,D,~]=eig(A);
        %find max eigenvalue
        maxEig=max(abs(diag(D)));
        if isreal(maxEig)
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
            A=X*DD/X;
            A=real(A);      %Extract real parts to eliminate roudnoff errors
        end
        %Check Anew has correct eigenvalues
        Dnew=sort(eig(A));
        if isreal(Dnew(end)) || ~isreal(Dnew(end)*Dnew(end-1))
            error("A does not have imaginary maximum eigenvalues")
        end
    case 'symmetric'
        A=rand(size);
        A=A'*A;
    case 'SPD'
        A=rand(size);A=A'*A;
        A=A+size*eye(size);
    case 'general'
        A=rand(size);
    case 'DiagDom'
        A=rand(size)+size*eye(size);
    case 'SymTriDiag'
        d=rand(size,1);
        dl=rand(size-1,1);
        A=diag(dl,-1)+diag(d,0)+diag(dl,1);
    otherwise
        error('Matrix Type not recognized')
end
end