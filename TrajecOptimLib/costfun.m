function J = costfun(X,CONSTANTS)

    X = matrify(X,CONSTANTS);
    J = 0;

    % dynamics
    J_dyn = zeros(CONSTANTS.Nv,1);
    for i = 1:CONSTANTS.Nv
        J_dyn(i) = CONSTANTS.costfun_single(X(:,:,i),CONSTANTS);
    end    
    J = J + sum(J_dyn);

    % obstacles
    if ~isempty(CONSTANTS.obstacles) && CONSTANTS.uselogbar
        J_obs = cell(CONSTANTS.Nv,1);
        for i = 1:CONSTANTS.Nv
            J_obs{i} = arrayfun(@(z)logbarrierfunc(0.1,z,CONSTANTS.usesigma), -env_collision_avoidance(X(:,:,i),CONSTANTS));
        end
        J = J+sum(cell2mat(J_obs));
    end
    if ~isempty(CONSTANTS.obstacles_circles) && CONSTANTS.uselogbar
        J_obs = cell(CONSTANTS.Nv,1);
        for i = 1:CONSTANTS.Nv
            J_obs{i} = arrayfun(@(z)logbarrierfunc(0.1,z,CONSTANTS.usesigma), -env_collision_avoidance_circles(X(:,:,i),CONSTANTS));
        end
        J = J+sum(cell2mat(J_obs));

    end
    
    % inter vehicles
    if CONSTANTS.Nv>1  && true
        J_barrier = cell(nchoosek(CONSTANTS.Nv,2),1);
        it = 1;
        for i = 1:CONSTANTS.Nv
            for j = i+1:CONSTANTS.Nv
                J_barrier{it} = arrayfun(@(z)logbarrierfunc(0.5,z,CONSTANTS.usesigma),-vehicle_collision_avoidance_isaac(X(:,1:2,i),X(:,1:2,j),CONSTANTS));          
                it = it+1;
            end
        end
        J = J+sum(cell2mat(J_barrier));
    end

end