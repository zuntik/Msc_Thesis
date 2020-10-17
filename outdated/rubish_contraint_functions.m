
function c = env_collision_avoindance_more_or_less_wrong(X,CONSTANTS)

    c = [-1;-1 ];
    
    if isempty(CONSTANTS.obstacles)
        return
    end
    
    N = CONSTANTS.N;
    
    cpts = X(:,1:2);
    curvepoints = CONSTANTS.EvalMat*X(:,1:2);

    for it = 1:size(CONSTANTS.obstacles,3)
        obstacle = CONSTANTS.obstacles(:,:,it);
        p1 = ConvexPolygon(obstacle);    
        upmin = -1;
        sidemin = -1;

        allpoligons = cell(N);
        for i = 1:N
            if i == 1
                vertices = [cpts(i:i+1,:);curvepoints(i+1,:)];
            elseif i == N
                vertices = [cpts(i:i+1,:);curvepoints(i,:)];
            else
                vertices = [cpts(i:i+1,:);curvepoints(i:i+1,:)];
            end
            vertices = uniquetol(vertices,'ByRows',true);
            allpoligons{i} = sortvertices(vertices);
            p2 = ConvexPolygon(vertices);
            [intersection,simplex] = gjk2(p1,p2);
            if intersection
                support=@(d) p1.support(d)-p2.support(-d);
                vec_up = directed_epa(support, simplex, [0;1]);
                vec_side = directed_epa(support, simplex, [1;0]);
                if vec_up(2) > upmin 
                    upmin = vec_up(2);
                end
                if vec_side(1) > sidemin
                    sidemin = vec_side(1);
                end
            end
        end
    end
    
    c = [upmin; sidemin];
    
end

function c = env_collision_avoindance_fair_wrong(X,CONSTANTS)

    c = [-1;-1 ];
    
    if isempty(CONSTANTS.obstacles)
        return
    end
    
    N = CONSTANTS.N;
    
    cpts = X(:,1:2);
    curvepoints = CONSTANTS.EvalMat*X(:,1:2);

    for it = 1:size(CONSTANTS.obstacles,3)
        obstacle = CONSTANTS.obstacles(:,:,it);
        p1 = ConvexPolygon(obstacle);    
        allpoligons = cell(N,1);
        intersecting_indexes = false(N,1);
        for i = 1:N
            if i == 1
                vertices = [cpts(i:i+1,:);curvepoints(i+1,:)];
            elseif i == N
                vertices = [cpts(i:i+1,:);curvepoints(i,:)];
            else
                vertices = [cpts(i:i+1,:);curvepoints(i:i+1,:)];
            end
            vertices = uniquetol(vertices,'ByRows',true);
            allpoligons{i} = vertices;
            p2 = ConvexPolygon(vertices);
            intersecting_indexes(i) = gjk2(p1,p2);
        end
        if ~any(intersecting_indexes)
            continue
        end
        vertices = cell2mat(allpoligons(intersecting_indexes));
        p2 = ConvexPolygon(uniquetol(vertices,'ByRows',true));
        [~,simplex] = gjk2(p1,p2);
        support=@(d) p1.support(d)-p2.support(-d);
        vec_up = directed_epa(support, simplex, [0;1]);
        vec_side = directed_epa(support, simplex, [1;0]);
    end
    
    if ~any(intersecting_indexes)
        c = [-1;-1];
    else
        c = [vec_up(2); vec_side(1)];
    end
    
end

function c = env_collision_avoindance_wrong(X,CONSTANTS)
    
    if isempty(CONSTANTS.obstacles)
        c = [];
        return
    end
    
    N = CONSTANTS.N;
    
    cpts = X(:,1:2);
    curvepoints = CONSTANTS.EvalMat*X(:,1:2);

    c = -1.*ones(N,size(CONSTANTS.obstacles,3));
    for it = 1:size(CONSTANTS.obstacles,3)
        obstacle = CONSTANTS.obstacles(:,:,it);
        pobs = ConvexPolygon(obstacle);
        intersecting_indexes = recursive_collision_test(cpts,curvepoints,1,N+1,pobs);
        if ~any(intersecting_indexes)
            continue
        end
        for i = 1:N
            if intersecting_indexes(i)&&intersecting_indexes(i+1)
                pcurv = ConvexPolygon([ cpts(i:i+1,:); curvepoints(i:i+1,:)]);
                [~,simplex] = gjk2(pobs,pcurv);
                support=@(d) pobs.support(d)-pcurv.support(-d);
                vec_up = directed_epa(support,simplex,[0;1]);
                c(i,it) = vec_up(2);
            end
        end
    end
    c = c(:);
    
end

function c = env_collision_avoindance_good_wrong(X,CONSTANTS)

    c = [-1;-1 ];
    
    if isempty(CONSTANTS.obstacles)
        return
    end
    
    N = CONSTANTS.N;
    
    cpts = X(:,1:2);
    curvepoints = CONSTANTS.EvalMat*X(:,1:2);

    for it = 1:size(CONSTANTS.obstacles,3)
        obstacle = CONSTANTS.obstacles(:,:,it);
        pobs = ConvexPolygon(obstacle);    
        intersecting_indexes = recursive_collision_test(cpts,curvepoints,1,N,pobs);
        if ~any(intersecting_indexes)
            continue
        end
        vertices = [ cpts(intersecting_indexes,:); curvepoints(intersecting_indexes,:)];
        pcurv = ConvexPolygon(vertices);
        [~,simplex] = gjk2(pobs,pcurv);
        support=@(d) pobs.support(d)-pcurv.support(-d);
        vec_up = directed_epa(support, simplex, [0;1]);
        %vec_side = directed_epa(support, simplex, [1;0]);
    end
    
    if ~any(intersecting_indexes)
        %c = [-1;-1];
        c  = -1;
    else
        %c = [vec_up(2); vec_side(1)];
        c = vec_up(2);
    end
    
end

function indexes = recursive_collision_test(ctrl_pts,curv_pts,l,r,p_obs)

    p_curv = ConvexPolygon([ctrl_pts(l:r,:);curv_pts(l:r,:)]);
    indexes = false(r-l+1,1);
    if ~gjk2(p_obs,p_curv)
        return
    else
        if r-l == 1
            indexes = true(2,1);
            return
        else
            halfsize = floor((r-l)/2)+1;
            midp = l+halfsize-1;% middle point
            indexes(1:halfsize) = indexes(1:halfsize)|recursive_collision_test(ctrl_pts,curv_pts,l,midp,p_obs);
            indexes(halfsize:end) = indexes(halfsize:end)|recursive_collision_test(ctrl_pts,curv_pts,midp,r,p_obs);
        end
    end

end
