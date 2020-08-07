function circs = surroundhull(x)
    % x are 2D obstacles
    circs = zeros(size(x,3),3);
    for i = 1:size(x,3)
        cx = convhull(x(:,:,i),'Simplify',true);
        [centrx,centry] = centroid(polyshape(x(cx(1:end-1),:,i)));
        r = max(sqrt( sum((x(cx(1:end-1),:)-[centrx,centry]).^2,2)));
        circs(:,i) = [ centrx, centry, r];
    end
end