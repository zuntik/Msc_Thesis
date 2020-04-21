function [new_p] = BernsteinSum(p1,p2)

if size(p1,1)==1 && size(p1,2) ~= 1
    p1 = p1';
    p2 = p2';
end

[ m, dim1 ] = size(p1);
[ n, dim2 ] = size(p2);

if dim1 ~= dim2
    error('dimentions don''t match');
end

if m>n
    p2 = BernsteinDegrElev(p2,m-1);
else
    p1 = BernsteinDegrElev(p1,n-1);
end

new_p = p1+p2;

%     if size(p1,1)==1 && size(p1,2) ~= 1
%         p1 = p1';
%         p2 = p2';
%     end
% 
%     if size(p1,1)<size(p2,1)
%         tmp = p1;
%         p1 = p2;
%         p2 = tmp;
%     end
%     
%     m = size(p1,1)-1;
%     n = size(p2,1)-1;
%     
%     new_p = zeros(m+1,size(p1,2));
%     
%     for i = 0:m
%         new_p(i+1,:) = p1(i+1,:);
%         for j = max([0,i-m+n]):min(n,i)
%             new_p(i+1,:) = new_p(i+1,:) + (nchoosek(n,j)*nchoosek(m-n,i-j)/nchoosek(m,i)).*p2(j+1,:);
%         end
%     end    
    
end