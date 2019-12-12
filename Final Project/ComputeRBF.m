function [fvalue] = ComputeRBF(basis,Beta,x,xi)
r=sqrt(sum((x-xi).^2));
fvalue=basis(r)'*beta;
end

