function [c,ceq] = nonlcon(X,CONSTANTS)

    X = matrify(X(:),CONSTANTS);

    c = [];
    ceq = [];

    % dynamics
    ceq_dyn = cell(CONSTANTS.Nv,1);
    c_dyn = cell(CONSTANTS.Nv,1);
    for i = 1:CONSTANTS.Nv
        [c_i,ceq_i] = CONSTANTS.dynamics(X(:,:,i),CONSTANTS);
        ceq_dyn{i} = ceq_i;
        c_dyn{i} = c_i;
    end
    ceq = [ceq; cell2mat(ceq_dyn)];
    c = [c; cell2mat(c_dyn)];
    
    % obstacles
    if ~isempty(CONSTANTS.obstacles)% && ~CONSTANTS.uselogbar
        c_env_avoidance = cell(CONSTANTS.Nv,1);
        for i = 1:CONSTANTS.Nv
            c_env_avoidance{i} = env_collision_avoidance(X(:,:,i),CONSTANTS);
        end    
        c = [c; cell2mat(c_env_avoidance)];
    end
    if ~isempty(CONSTANTS.obstacles_circles)% && ~CONSTANTS.uselogbar
        c_env_avoidance = cell(CONSTANTS.Nv,1);
        for i = 1:CONSTANTS.Nv
            c_env_avoidance{i} = env_collision_avoidance_circles(X(:,:,i),CONSTANTS);
        end
        c = [c; cell2mat(c_env_avoidance)];
    end

    % inter vehicles
    if CONSTANTS.Nv>1% && ~CONSTANTS.uselogbar
        c_vehicle_avoidance = cell(nchoosek(CONSTANTS.Nv,2),1);
        it = 1;
        for i = 1:CONSTANTS.Nv
            for j = i+1:CONSTANTS.Nv
                %c_vehicle_avoidance{it} = vehicle_collision_avoidance_bad(X(:,1:2,i),X(:,1:2,j),CONSTANTS);
                c_vehicle_avoidance{it} = vehicle_collision_avoidance_isaac(X(:,1:2,i),X(:,1:2,j),CONSTANTS);
                it = it+1;
            end
        end
        c = [c;cell2mat(c_vehicle_avoidance)];
    end

end