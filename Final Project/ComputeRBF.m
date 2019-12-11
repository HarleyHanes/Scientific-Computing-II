function [fvalue] = ComputeRBF(basis,Beta,x,xi)
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here
%Compute Vector of Radi

r=sqrt(sum((x-xi).^2));
fvalue=basis(r)'*beta;
end

