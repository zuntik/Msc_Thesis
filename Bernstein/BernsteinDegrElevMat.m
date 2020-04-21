function [mat] = BernsteinDegrElevMat(n,m)
%BERNSTEINDEGRELEVMAT Summary of this function goes here
%   Detailed explanation goes here

mat = zeros(m+1,n+1);

for i = 0:m
    for j = [ max(0, i-m+n) : min(n,i) ]
        %mat(i+1,j+1) = nchoosek(n,j)*nchoosek(m-n,i-j)/nchoosek(m,i);
        mat(i+1,j+1) = nchoosek(i,j)*nchoosek(m-i,n-j);
    end
end
mat = mat./nchoosek_mod(m,n);

end