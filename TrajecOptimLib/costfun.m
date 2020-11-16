function J = costfun(X,constants)

    J = 0;
    if constants.uselogbar
        c = nonlcon_ineq(X, constants);
        J = J + sum(logbarrierfunc(0.1, -c, constants.usesigma));
        if constants.useeqlogbar
            ceq = nonlcon_eq(X, constants);
            %J=J+sum(logbarrierfunc(0.01,-abs(ceq),constants.usesigma));
            %J=J+sum(logbarrierfunc(0.01,1e-5-ceq.^2,constants.usesigma));
            %J=J+logbarrierfunc(0.1,-sum(1e-5-ceq.^2),constants.usesigma);
            %J = J + sum(ceq.^2);
            J=J+sum(...
                logbarrierfunc(0.1,-ceq,constants.usesigma)+...
                logbarrierfunc(0.1,ceq,constants.usesigma)...
                );
        end
        %J = J*10e-3;
    end
    X = matrify(X,constants);
    % dynamics
    J_dyn = zeros(constants.Nv,1);
    for i = 1:constants.Nv
        J_dyn(i) = constants.costfun_single(X(:,:,i),constants);
    end    
    J = J + sum(J_dyn);


%     % obstacles
%     if ~isempty(constants.obstacles) && constants.uselogbar
%         J_obs = cell(constants.Nv,1);
%         for i = 1:constants.Nv
%             J_obs{i} = arrayfun(@(z)logbarrierfunc(0.1,z,constants.usesigma), -env_collision_avoidance(X(:,:,i),constants));
%         end
%         J = J+sum(cell2mat(J_obs));
%     end
%     if ~isempty(constants.obstacles_circles) && constants.uselogbar
%         J_obs = cell(constants.Nv,1);
%         for i = 1:constants.Nv
%             J_obs{i} = arrayfun(@(z)logbarrierfunc(0.1,z,constants.usesigma), -env_collision_avoidance_circles(X(:,:,i),constants));
%         end
%         J = J+sum(cell2mat(J_obs));
% 
%     end
%     
%     % inter vehicles
%     if constants.Nv>1  && true
%         J_barrier = cell(nchoosek(constants.Nv,2),1);
%         it = 1;
%         for i = 1:constants.Nv
%             for j = i+1:constants.Nv
%                 J_barrier{it} = arrayfun(@(z)logbarrierfunc(0.5,z,constants.usesigma),-vehicle_collision_avoidance_isaac(X(:,1:2,i),X(:,1:2,j),constants));          
%                 it = it+1;
%             end
%         end
%         J = J+sum(cell2mat(J_barrier));
%     end

end