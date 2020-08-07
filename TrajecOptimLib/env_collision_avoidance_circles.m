function c = env_collision_avoidance_circles(X,CONSTANTS)
    c = [];
    if isempty(CONSTANTS.obstacles_circles)
        return
    end
    c = zeros(1,size(CONSTANTS.obstacles_circles,1));
    for i = length(c)
        c(i) = CONSTANTS.obstacles_circles(i,3) - min(sqrt(sum((CONSTANTS.BigElevMat*(X(:,1:2)-CONSTANTS.obstacles_circles(i,1:2))).^2,2)));
    end
end