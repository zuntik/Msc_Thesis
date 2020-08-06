function d = MinTimeDistanceBernstein(X,varargin)

    % --- Input parsing
    p = inputParser;

    addParameter(p, 'tolerance', 1e-5, @isnumeric)

    parse(p, varargin{:})
    tolerance = p.Results.tolerance;
    % --- End input parsing

    % X is size [N+1,dim,Ncuves]

    [N,dim,Ncurves] = size(X);
    if Ncurves<2
        error('no point in running this function with just one curve')
    end

    if dim~=2
        error('this function only supports 2D');
    end

    N = N-1;

    maxiter = 1000;

    % d is a maximum value for min dist
    d = min([   closest_points(reshape(X(1,  :,:),[dim,Ncurves]).'),...
                closest_points(reshape(X(end,:,:),[dim,Ncurves]).')]);

    X = {X};

    X = splitcurves(X);

    % use increaseingly sized list instead of recusion
    for i = 1:maxiter
%         plotstuff(X);
        d = min([testlastpoints(X) d]);
        [X, maxrad] = activecurves(X,d);
        if maxrad < tolerance
            break
        end
        X = splitcurves(X);

    end
    if i == maxiter
        warning('max interations met');
    end

end

function [c,pca] = surroundhull(x)

    cx = convhull(x,'Simplify',true);
    pca = polyshape(x(cx(1:end-1),:));
    [centrx,centry] = centroid(pca);
    r = max(sqrt( sum((x(cx(1:end-1),:)-[centrx,centry]).^2,2)));
    c = [ centrx, centry, r];
    
end

function [Xout,maxrad] = activecurves(X,d)
    
    maxrad = 0;
    Xout = cell(length(X),1);
    for i = 1:length(X)
        Xcur = X{i};
        circles = zeros(size(Xcur,3),3);
        for j = 1:size(Xcur,3)
            circles(j,:)=surroundhull(Xcur(:,:,j));
        end
        intersect_indexes = circles_intersect(circles + [ 0 0 d/2]);
        %intersect_indexes = polygons_intersect(Xcur,d);
        maxrad = max([maxrad; circles(:,3)]);
        Xout{i} = Xcur(:,:,intersect_indexes);
    end
    Xout = Xout(~cellfun(@isempty,Xout));
    
end

function intersect_indexes = polygons_intersect(X,d)
    % peform 
end

function X = splitcurves(Xin)
    % spit the cuves in half the time.
    Ns = length(Xin); % number of segments
    X = cell(Ns*2,1);
    for i = 1:Ns
        x = Xin{i};
        for j = 1:size(x,3)
            x_cut = my_deCasteljau(x(:,:,j),0.5);
            X{i}(:,:,j) = x_cut(:,:,1);
            X{i+Ns}(:,:,j) = x_cut(:,:,2);
        end
    end
    
end

function d = testlastpoints(X)
    % testes the last points of the first half of "cubes"

    d = Inf;
	for i = 1:floor(length(X)/2)
        x = X{i};
        endpoints = reshape(x(end,:,:),[size(x,2),size(x,3)]).';
		d = min([closest_points(endpoints), d]);
	end
end

function plotstuff(X)
    figure, axis equal, hold on
    for i = 1:length(X)
        xcur = X{i};
        for j = 1:size(xcur,3)
            [~,pca] = surroundhull(xcur(:,:,j));
            plot(pca);
        end
    end
end