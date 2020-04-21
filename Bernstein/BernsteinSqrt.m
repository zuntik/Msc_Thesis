function [ p1, p2 ] = BernsteinSqrt(s)

    if size(s,1)==1 && size(s,2) ~= 1
        s = s';
    end

    n2 = size(s,1) - 1;

    if rem(n2,2) ~= 0
        n2 = n2+1;
        s = BernsteinDegrElev(s,n2);
    end

    n = n2/2;

    p = zeros(n+1,size(s,2),2);
    p(1,:,1) =  sqrt(s(1,:));
    p(1,:,2) = -sqrt(s(1,:));

    for i = 1:n
        denom = 2 .* nchoosek(n,i) .* p(1,:,:);
        for j = 1:i-1
            denom = denom + nchoosek(n,j) .* nchoosek(n,i-j) .* p(j+1,:,:) .* p(i-j+1,:,:);
        end
        p(i+1,:,:) = nchoosek(n2,i) .* s(i+1,:) ./ denom;
    end

    p1 = p(:,:,1);
    p2 = p(:,:,2);

end