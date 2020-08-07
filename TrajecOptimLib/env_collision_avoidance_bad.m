function c = env_collision_avoidance_bad(X,CONSTANTS)

    c = zeros(2,size(CONSTANTS.obstacles,3));
    for i = 1:size(CONSTANTS.obstacles,3)
        obstacle = CONSTANTS.obstacles(:,:,i);        
        [dist, t, pt] = BernsteinMinDist2Shape(X(:,1:2),ConvexPolygon(obstacle));
        if dist>0 || length(t) == 1
            c(:,i) = [ -1; -1];
        else
            if length(t) == 2 && false
                cptsleftcut = my_deCasteljau(X(:,1:2),t(1));
                cptsrightcut = my_deCasteljau(cptsleftcut(:,:,2),t(2));
                cpts_inside = cptsrightcut(:,:,1);
            else
                cpts_inside = pt;
            end
            p2 = ConvexPolygon(cpts_inside);
            p1 = ConvexPolygon(obstacle);
            [~,~,simplex] = gjk2D(p1,p2,'collcheckonly',true);
            support=@(d) p1.support(d)-p2.support(-d);
            vec_up = directed_epa(support, simplex, [0;1]);
            vec_side = directed_epa(support, simplex, [1;0]);
            c(:,i) = [vec_up(2); vec_side(1)];
        end
    end
    c = c(:);
    
end
