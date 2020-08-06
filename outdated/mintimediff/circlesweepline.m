function indexes = circlesweepline(circles)
% circles is mat with 3 columns: x,y,r (and z in future)

	Nc = size(circles,1); % number of circles

	cp = cell(2*Nc,1); % circle points

	for i = 1:Nc
		cp{i*2-1} = struct('id', i, 'isLeft', true,  'x', circles(i,1)-circles(i,3));
		cp{i*2}   = struct('id', i, 'isLeft', false, 'x', circles(i,1)+circles(i,3));
	end

	x_min = circles(:,1) - circles(:,3);
	x_max = circles(:,1) + circles(:,3);

	[sorted_x,indexes] = sort([x_min;x_max]);

	% no point in testing pairs that have been tested before
	tested_combinations = rot90(diag(ones(Nc,1)))==1;
	indexes = false(Nc,1);

	t = cell(); %ordered list

	for i = 1:2*Nc
		if cp{i}.isLeft
			listinsert(t,cp{i}.down);
			listinsert(t,cp{i}.up);
			prev = listfindprev(t,cp{i});
			if intersects(cp{i},prev)
				% paint cur and prev
				indexes(cp{i}.id) = true;
				indexes(prev.id) = true;
				% tag matrix as listed
				
			end
			next = listfindnext(t,cp{i});
			if intersects(cp{i},next)
				% paint cur and next
				indexes(cp{i}.id) = true;
				indexes(next.id) = true;
			end
		else
			if intersects(listfindprev(t,cp{i}),listfindnext(t,cp{i}))
				% paint them
			end
			t.delete() % delete both bot and top
		end
	end


end
