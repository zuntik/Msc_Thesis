function c = env_collision_avoidance_circles_venanzio(X,constants)
    persistent nchoosek_mat N
    if isempty(nchoosek_mat) || constants.N ~= N
        N = constants.N;
        nchoosek_mat = nchoosek_mod_mat(2*N+3);
    end
    c = [];
    if isempty(constants.obstacles_circles)
        return
    end
    c = cell(1,size(constants.obstacles_circles,1));
    for i = length(c)
        cir = constants.obstacles_circles(i,:);
        dist2obs2 = BernsteinPow(X(:,1)-cir(1),2, nchoosek_mat) + BernsteinPow(X(:,2)-cir(2),2, nchoosek_mat);
        c{i} = -dist2obs2 + cir(3)^2;        
    end
    c = cell2mat(c);
end