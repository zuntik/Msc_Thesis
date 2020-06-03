function indexes = grahamscan(A)
% A is [N x 2]  array
% grahamscan - perform the Graham Scan to find the minimal convex hull
% of a polygon. 
% It will also be able to "sort" the polygons
% Assumes that A has no repeated points, can be achieved with 

%% First remove duplicate points from A 

[A,ia,ic] = unique(A,'rows','stable');


%% Find the bottommost point 
ymin = A(1,2); argymin = 1;
for i = 2:size(A,1)
    y = A(i,2);
    % Pick the bottom-most or chose the left 
    % most point in case of tie 
    if y < ymin || ( ymin == y && A(i,1) < A(argymin,1) )
        ymin = A(i,2); argymin = i;
    end
end

%% Sort n-1 points with respect to the first point. 
% A point p1 comes before p2 in sorted output if p2  has larger polar
% angle (in counterclockwise direction) than p1 

bottomleftmost = A(argymin,:);
% A(argymin,:) = [];

downtoA = A - bottomleftmost;

norms_squared = sqrt(downtoA(:,1).^2+downtoA(:,2).^2);
allcosines = downtoA(:,1)./norms_squared;
allcosines(argymin)=-2;
% reinclude the most bottomleft point where it was 
% allcosines = [allcosines(1:argymin-1); -2;  allcosines(argymin:end) ];
% A = [ A(1:argymin-1,:); bottomleftmost; A(argymin:end,:)];
% downtoA = [downtoA(1:argymin-1,:); [0 0]; downtoA(argymin:end,:)];

[uniquecosines,~,iccosines] = unique(allcosines);

candidates = zeros(size(A,1),1);

for i = 1:length(uniquecosines)
    indexes = find(iccosines==i);
    [~,argmaxnorm] = max(norms_squared(indexes));
    candidates(indexes(argmaxnorm))=1;
end

% now sort!
% notice the bottom left point will most definitely be first on the list
[~,sortedpos] = sort(uniquecosines);

A =A(find(candidates),:);

if length(uniquecosines) <3
    indexes=ones(length(uniquecosines),1);
    indexes = recoverindexes(indexes,candidates,ic,ia);
    return
    
end

%% perform scan
% if oritentation > 0 clockwise, elseif 0 colinear else anticlockwise
orientation = @(p1,p2,p3) (p2(2)-p1(2))*(p3(1)-p2(1))-(p2(1)-p1(1))*(p3(2)-p2(2));

% create the stack that contains the indexes of points in hull
l = java.util.LinkedList();
l.push(sortedpos(1));
l.push(sortedpos(2));
l.push(sortedpos(3));

for i = 4:size(A,1)

    while true
        p1 = A(l.get(1),:);
        p2 = A(l.get(0),:);
        p3 = A(sortedpos(i),:);
        ori = orientation(p1,p2,p3);
        if ori<=0
            l.pop();
        else
            break
        end
    end

    l.push(sortedpos(i));

end
%% finish up by creating the vector of indexes
% ... of points belonging to the convex hull

indexes = zeros(size(A,1),1); % is part of the convex hull

it = l.listIterator();
while it.hasNext
    indexes(it.next())=1;
end

indexes = recoverindexes(indexes,candidates,ic,ia);



function indexes = recoverindexes(indexes,candidates,ic,ia)
    % bring back the points with same cosine
    id = indexes;
    indexes = zeros(nnz(candidates),1);
    indexes(find(candidates)) = id;

    % bring back the repeated points
    id = indexes;
    % assume that duplicated points do not belong in the convex hull
    indexes = zeros(length(ic),1);
    indexes(ia) = id;
end

end