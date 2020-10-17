function c = env_collision_avoidance_circles(X,constants)
    c = [];
    if isempty(constants.obstacles_circles)
        return
    end
    c = zeros(1,size(constants.obstacles_circles,1));
    for i = length(c)
        c(i) = constants.obstacles_circles(i,3) - min(sqrt(sum((constants.BigElevMat*(X(:,1:2)-constants.obstacles_circles(i,1:2))).^2,2)));
    end
end