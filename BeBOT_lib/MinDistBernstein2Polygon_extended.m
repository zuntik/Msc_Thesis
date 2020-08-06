function [dist, t, pt] = MinDistBernstein2Polygon_extended(cpts, poly, varargin)
%BEZ2POLY Finds the minimum distance between a Bezier curve and a polytope
% or returns where curve intersects
%   INPUTS
%   cpts: Dx(N+1) Dth dimensional array of N+1 control points defining the
%       first Nth order curve.
%   poly: DxM Dth dimensional array of M vertices defining a polygon.
%   epsilon: Optional parameter, desired tolerance to meet. Default is 1e-6
%
%   OUTPUTS
%   dist: Minimum distance between the two curves
%   t: t value at which the minimum distance happens on the curve.
%       Note that t will be on the range [0, 1]
%   pt: Point at which the minimum distance happens on the polytope

% Written by Calvin Kielas-Jensen
% Exdended by Thomas Berry

[m,n] = size(cpts);
if m == 2
    cpts = [cpts;zeros(1,n)];
end
[m,n] = size(poly);
if m == 2
    poly = [poly;zeros(1,n)];
end

% --- Input parsing
p = inputParser;

addParameter(p, 'epsilon', 1e-6, @isnumeric)

parse(p, varargin{:})
epsilon = p.Results.epsilon;
% --- End input parsing

alpha = inf;
deg = size(cpts, 2) - 1;

[dist, t, pt, intersectpoints] = bern2polyR(cpts, poly, alpha, deg, epsilon, 0, 1, 0, false);

if ~isempty(intersectpoints)
    t = intersectpoints(:).';
    pt = BernsteinPolyStable(cpts,t(:).');
    dist = -1;
end

if m == 2
    pt = pt(1:2,:);
end

end



function [alpha, t, pt, intersectpoints] = bern2polyR(...
    cpts, poly, alpha, deg, epsilon, tL, tH, count,ispartiallyin)
%bern2polyR Recursive call for bern2polyR

% Maximum iteration check
count = count + 1;
if count > 1000
    t = 0.5*(tH + tL);
    intersectpoints=[]; pt=[]; alpha=[];
    warning('B2PR: Maximum number of iterations met.')
    return
end

% if all of the control points are in the polygon, no min and no intersect
% if iscontained(poly.',cpts.')
%     intersectpoints=[]; pt = []; t = []; alpha = [];
%     return
% end
[cpts_in_obs,cpts_edge_obs] = inpolygon(cpts(1,:),cpts(2,:),poly(1,:),poly(2,:));
if all(cpts_in_obs) && ~any(cpts_edge_obs)
    intersectpoints=[]; pt = []; t = []; alpha = [];
    return
end

[ub, tlocal, pt, partiallyin] = upper_bound_poly(cpts, poly);

% if the curve is partially in the object no point in min dist
if ~partiallyin
    [lb, ~, ~, simplex] = gjk(cpts, poly);
else
    lb = 0;
end

% if we already know that the curve intersects at least once then no need
% to bother to calculate the minimum distance
if ispartiallyin && lb > 0
    intersectpoints=[]; pt = []; t = []; alpha = [];
    return
end

% if partiallyin then ub becomes the distence between tL and tH
if partiallyin && ub < epsilon
    intersectpoints = tL;
    %pt = cpts(:,1);
    alpha = [];
    t = [];
    return
end

% if lb >= alpha*(1-epsilon) % Relative tolerance
if ~partiallyin && (lb >= alpha - epsilon) % Absolute tolerance
    intersectpoints=[];
    t = (1-tlocal)*tL + tlocal*tH;
    % pt has already been determined by the upper bound function
    return
end

% if partially in then the values of alpha and ub mean nothing
if ~partiallyin && ub < alpha
    alpha = ub;
    t = (1-tlocal)*tL + tlocal*tH;
else
    t = [];
end

% now determine where to split
if ~partiallyin && lb > 0
    tdiv = calculatesplitpoint(cpts,simplex,deg);
else
    % If there was a collision, divide at 0.5
    tdiv = 0.5; 
end

% Split the curve
cptsP = deCasteljau(cpts, tdiv);
cptsPA = cptsP(:, 1:deg+1);
cptsPB = cptsP(:, deg+1:end);

tlen = tH - tL;

% Recursively call and keep the minimum value along with its
% corresponding t values
% the followining 
intersectpoints = [];  % this isn't python but this array is meant to be though of as a list of times

[newAlpha, newT, newPt, intersectpointsA] = bern2polyR(cptsPA, poly, alpha, deg, ...
                    epsilon, tL, tL + tdiv*tlen, count, partiallyin);

if ~isempty(newAlpha) && newAlpha < alpha
    alpha = newAlpha;
    t = newT;
    pt = newPt;
    intersectpoints=[];
end
    
if ~isempty(intersectpointsA)
    intersectpoints = [intersectpoints; intersectpointsA];
    t = [];
end

[newAlpha, newT, newPt, intersectpointsB] = bern2polyR(cptsPB, poly, alpha, deg, ...
                    epsilon, tL + tdiv*tlen, tH, count, partiallyin);

if ~isempty(newAlpha) && newAlpha < alpha
    alpha = newAlpha;
    t = newT;
    pt = newPt;
    intersectpoints=[];
end
    
if ~isempty(intersectpointsB)
    intersectpoints = [intersectpoints; intersectpointsB];
    t = [];
end

end



function tdiv = calculatesplitpoint(cpts,simplex,deg)

% Convert the struct to a cell array
vertices = cell(1, simplex.len);
baryVals = zeros(1, simplex.len);
if simplex.len >= 1
    vertices{1} = simplex.Apts;
    baryVals(1) = simplex.Au;
end
if simplex.len >= 2
    vertices{2} = simplex.Bpts;
    baryVals(2) = simplex.Bu;
end
if simplex.len >= 3
    vertices{3} = simplex.Cpts;
    baryVals(3) = simplex.Cu;
end
if simplex.len >= 4
    vertices{4} = simplex.Dpts;
    % There should not be a Du value since that would be for a 4D shape
end

% Find the index of where the simplex vertices are in cpts
approx = @(x, y) abs(x-y) < 1000000*eps; 
idx = zeros(1, simplex.len);
for i = 1:simplex.len
    idx(i) = find(all(bsxfun(approx, cpts, vertices{i}(:, 1)), 1), 1)-1;
end

% Find the point at which to split the curves
if sum(baryVals) > 1
    tdiv = 0.5;
else
    % By finding a good place to split, we significantly speed up the
    % convergence time
    tdiv = sum(idx.*baryVals / deg);
    if (tdiv == 0) || (tdiv == 1)
        tdiv = 0.5;
    end
end

end



function [ub, t, pt, partiallyin] = upper_bound_poly(cpts, poly)
%UPPER_BOUND_POLY Finds the upper bound on the distance between a Bezier
%   curve and a polytope.
%
%   Detailed explanation goes here

[d1, ~, p1, simplex1] = gjk(cpts(:,   1), poly);
[d2, ~, p2, simplex2] = gjk(cpts(:, end), poly);

partiallyin = simplex1.collision || simplex2.collision;
if partiallyin
    ub = sum((cpts(:,1)-cpts(:,2)).^2);
    t = 0; pt = cpts(:,1);
else
    if d1 < d2
        t = 0;
        pt = p1;
        ub = d1;
    else
        t = 1;
        pt = p2;
        ub = d2;
    end
end


end
