function indexes = circles_intersect(circles)
% circles is mat with 3 columns: x,y,r (and z in future)

	Nc = size(circles,1); % number of circles
	cp = cell(2*Nc,1); % circle points

	for i = 1:Nc
		cp{i} = struct('id', i, 'isLeft', true, 'c', circles(i,:) ); %,  'x', circles(i,1)-circles(i,3));
		cp{Nc+i}   = struct('id', i, 'isLeft', false, 'c', circles(i,:) ); % , 'x', circles(i,1)+circles(i,3));
	end

	x_min = circles(:,1) - circles(:,3);
	x_max = circles(:,1) + circles(:,3);

	[~,sorted_x] = sort([x_min;x_max]);

	indexes = false(Nc,1);

	l = [];

	for i = 1:2*Nc
		curpoint = sorted_x(i);
		if cp{curpoint}.isLeft
			ids = findintersections(cp{curpoint}.id,l, circles);
            indexes(ids) = true;
			l = listinsert(l,cp{curpoint}.id);
		else
			l = listremove(l,cp{curpoint}.id);
		end
	end


end

function l = listinsert(l, id)
	l=[l;id];
end


function l = listremove(l, id)
	l(l==id) = [];
end


function ids = findintersections(cp, l, circles)

	indexes = false(length(l),1);
	for i = 1:length(l)
        x1=circles(cp,1);y1 = circles(cp,2);r1 = circles(cp,3);
        x2=circles(l(i),1);y2 = circles(l(i),2);r2 = circles(l(i),3);        
		indexes(i) = indexes(i) || ccintersect(x1,y1,x2,y2,r1,r2);
    end
    ids = l(indexes);
    if ~isempty(ids)
        ids = [ids;cp];
    end
%     disp(indexes)
% 	indexes = l(indexes);

end

function intersect = ccintersect(x1, y1, x2, y2, r1, r2)

    distSq = (x1 - x2) .* (x1 - x2) + (y1 - y2) .* (y1 - y2);  
    radSumSq = (r1 + r2) .* (r1 + r2);  
	intersect = distSq<=radSumSq;
end

