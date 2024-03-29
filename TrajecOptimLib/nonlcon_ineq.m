function c = nonlcon_ineq(X,constants)

    X = matrify(X(:),constants);

    c = [];

    % dynamics
    c_dyn = cell(constants.Nv,1);
    for i = 1:constants.Nv
        [c_i,~] = constants.dynamics(X(:,:,i),constants);
        c_dyn{i} = c_i;
    end
    c = [c; cell2mat(c_dyn)];
    
    % obstacles
    if ~isempty(constants.obstacles)% && ~constants.uselogbar
        c_env_avoidance = cell(constants.Nv,1);
        for i = 1:constants.Nv
            c_env_avoidance{i} = env_collision_avoidance(X(:,:,i),constants);
        end    
        c = [c; cell2mat(c_env_avoidance)];
    end
    if ~isempty(constants.obstacles_circles)% && ~constants.uselogbar
        c_env_avoidance = cell(constants.Nv,1);
        for i = 1:constants.Nv
            c_env_avoidance{i} = env_collision_avoidance_circles(X(:,:,i),constants);
        end
        c = [c; cell2mat(c_env_avoidance)];
    end

    % inter vehicles
    if constants.Nv>1% && ~constants.uselogbar
        c_vehicle_avoidance = cell(nchoosek(constants.Nv,2),1);
        it = 1;
        for i = 1:constants.Nv
            for j = i+1:constants.Nv
%                 c_vehicle_avoidance{it} = vehicle_collision_avoidance(X(:,1:2,i),X(:,1:2,j),constants);
                c_vehicle_avoidance{it} = vehicle_collision_avoidance_isaac(X(:,1:2,i),X(:,1:2,j),constants);
                %c_vehicle_avoidance{it} = vehicle_collision_avoidance_fast(X(:,1:2,i),X(:,1:2,j),constants);
                it = it+1;
            end
        end
        c = [c;cell2mat(c_vehicle_avoidance)];
    end

end