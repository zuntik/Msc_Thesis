function [min_val,indexes] = closest_points(P)
    [~,sort_x] = sort(P(:,1));
    P = P(sort_x,:);
    
    [min_val,indexes] = closestUtil(P);
    %[min_val, indexes] = bruteForce(P);
    indexes = sort_x(indexes);
    
end

function [min_val, indexes] = closestUtil(P)

    if size(P,1) <= 3
        [min_val,indexes] = bruteForce(P);
        return
    end
    
    mid = floor(size(P,1)/2);
    midPoint = P(mid,:);
    
    [dl,indexes_left] = closestUtil(P(1:mid,:));
    [dr,indexes_right] = closestUtil(P(mid+1:end,:));
    
    if dl<=dr
        indexes = indexes_left;
        d = dl;
    else
        indexes = indexes_right + mid;
        d = dr;
    end
    
    indexes_strip = abs(P(:,1)-midPoint(1)) < d;
    strip = P( indexes_strip,:);    
    [min_strip, indexes_in_strip] = stripClosest(strip,d);
    indexes_strip(indexes_strip==1) = indexes_in_strip;

	if d < min_strip
		min_val = d;
	else
		min_val = min_strip;
		indexes = find(indexes_strip);
	end
    
end

function [min_val, indexes] = stripClosest(strip, d)
    
    min_val = d^2;
    [~, sorted_strip] = sort(strip(:,2));
	strip = strip(sorted_strip,:);
    indexes = [];
    for i = 1:size(strip,1)
        j = i+1;
        while j <= size(strip,1) && (strip(j,2)-strip(i,2))^2 <min_val
            min_val = sum((strip(i,:)-strip(j,:)).^2);
			indexes = [i,j];
			j = j+1;
        end
    end
    if isempty(indexes)
        min_val = Inf;
    else
        min_val = sqrt(min_val);
    end
    id = sorted_strip(indexes);
    indexes= false(size(strip,1),1);
    indexes(id) = true;
    
end

function [min_val,indexes] = bruteForce(P)
    min_val = Inf;
    for i = 1:size(P,1)
        for j = i+1:size(P,1)
            if sum((P(i,:)-P(j,:)).^2) < min_val
                min_val = sum((P(i,:)-P(j,:)).^2);
                indexes = [ i;j];
            end
        end
    end
    min_val = sqrt(min_val);
end
