function c = env_collision_avoidance(X,constants)

    c = [];
    if isempty(constants.obstacles) 
        return
    end
    c = zeros(2,size(constants.obstacles,3));
    for i = 1:size(constants.obstacles,3)
        obstacle = constants.obstacles(:,:,i);
        [dist, t, pt] = MinDistBernstein2Polygon_extended(X(:,1:2).', obstacle.');
        if dist>0 || length(t) == 1
            c(:,i) = [ -1; -1];
        else
            if length(t) == 2 && false
                cptsleftcut = deCasteljau(X(:,1:2).',t(1));
                cptsrightcut = deCasteljau(cptsleftcut(:,N+1:end),t(2));
                cpts_inside = cptsrightcut(:,1:N+1);
            else
                cpts_inside = pt.';
            end
            p2 = ConvexPolygon(cpts_inside);
            p1 = ConvexPolygon(obstacle);
            [~,~,simplex] = gjk2D(p1,p2,'collcheckonly',true);
            support=@(d) p1.support(d)-p2.support(-d);
            vec_up = directed_epa(support, simplex, [0;1]);
            %vec_side = directed_epa(support, simplex, [1;0]);
            %c(:,i) = [vec_up(2); vec_side(1)];
            c(:,i) = [vec_up(2); -1];   
        end
    end
    c = c(:);
    
end