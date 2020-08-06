clear all; %#ok<CLALL>

addpath('..\Bernstein');
addpath('..\BeBOT_lib');

% load parameters

parameters

xIn = init_guess(CONSTANTS);

options = optimoptions(@fmincon,'Display','Iter','Algorithm','sqp','MaxFunctionEvaluations',300000,'StepTolerance',eps,'MaxIterations',Inf);
% options = optimoptions(@fmincon,'Display','Iter','Algorithm','sqp');

xOut_uc = zeros(CONSTANTS.N-1,CONSTANTS.numvars,CONSTANTS.Nv);
CONSTANTS2 = CONSTANTS;
CONSTANTS2.Nv=1;
CONSTANTS2.obstacles = [];
CONSTANTS2.obstacles_circles = [];
for i = 1:CONSTANTS.Nv
    CONSTANTS2.xi = CONSTANTS.xi(i,:);
    CONSTANTS2.xf = CONSTANTS.xf(i,:);
    xOut_uc(:,:,i) = fmincon(@(x)costfun(x,CONSTANTS2),xIn(:,:,i),[],[],[],[],[],[],@(x)nonlcon(x,CONSTANTS2),options);
end
if ~isempty(CONSTANTS.obstacles) || CONSTANTS.Nv>1 || ~isempty(CONSTANTS.obstacles_circles)
    tic
    xOut = fmincon(@(x)costfun(x,CONSTANTS),xOut_uc,[],[],[],[],[],[],@(x)nonlcon(x,CONSTANTS),options);
    toc
else
    xOut = xOut_uc;
end

%%

[c,ceq] = nonlcon(xOut,CONSTANTS);
disp((c<0).')
disp(norm(ceq))
disp(ceq.')

%%

% xOut = [CONSTANTS.xi; reshape(xOut,[],CONSTANTS.numvars); CONSTANTS.xf];
xOut = matrify(xOut,CONSTANTS);

%%
figure, axis equal, hold on
for i = 1:CONSTANTS.Nv
    BernsteinPlot(xOut(:,1:2,i),CONSTANTS.T);
    [~,xy] = recoverplot(xOut(:,:,i),CONSTANTS.T);
    plot(xy(:,1),xy(:,2));
end

if ~isempty(CONSTANTS.obstacles) && true
    for i = 1:size(CONSTANTS.obstacles,3)
        plot(polyshape(CONSTANTS.obstacles(:,:,i)))
    end
end

if ~isempty(CONSTANTS.obstacles_circles)
    for i = 1:size(CONSTANTS.obstacles_circles,3)
        t = 0:0.001:2*pi; 
        centrx = CONSTANTS.obstacles_circles(i,1);
        centry = CONSTANTS.obstacles_circles(i,2);
        r = CONSTANTS.obstacles_circles(i,3);
        plot(polyshape([(centrx+r*cos(0:0.01:2*pi));(centry+r*sin(0:0.01:2*pi))].'));
    end
end

if CONSTANTS.Nv == 2
    figure
    X_diff = xOut(:,1:2,1)-xOut(:,1:2,2);
    X_diff_2 = BernsteinPow(X_diff,2);
    Mag = sqrt(sum(X_diff_2,2));
    BernsteinPlot(Mag,CONSTANTS.T);
end

if ~isempty(CONSTANTS.obstacles_circles) && size(CONSTANTS.obstacles_circles,3) == 1 && CONSTANTS.Nv == 1
    figure, grid on
    BernsteinPlot(sum((xOut(:,1:2)-CONSTANTS.obstacles_circles(i,1:2)).^2,2),CONSTANTS.T);
end

%%
function [t,xy] = recoverplot(X,T)

    v = @(t) BernsteinEval(X(:,4),T,t);
    w = @(t) BernsteinEval(X(:,5),T,t);
    
    [t,xy] = ode45(@(t,xy)odefunc(t,xy,v,w), [0 T], X(1,1:3));

    function dydt = odefunc(t,y,v,w)

        dydt = zeros(3,1);
        dydt(1) = v(t)*cos(y(3));%x
        dydt(2) = v(t)*sin(y(3));%y
        dydt(3) = w(t);%psi
        
    end

end

function [c,ceq] = dynamics5vars(X,CONSTANTS)

    DiffMat = CONSTANTS.DiffMat;
    
    xp = X(:,1);
    yp = X(:,2);
    psi = X(:,3);
    v = X(:,4);
    w = X(:,5);

    ceq = [
        DiffMat*xp - v.*cos(psi)
        DiffMat*yp - v.*sin(psi)
        DiffMat*psi - w
    ];
    %c=[-v+0.8;v-1.2 ;psi-pi;-psi-pi];
    %c = [-v+0.2;psi-pi;-psi-pi]; % this one is venanzio's
    %c = [ -v; psi-pi-pi/3; -psi-pi-pi/3];
    c = -v + 0.2 ;
    
end

function circs = surroundhull(x)
    % x are 2D obstacles
    circs = zeros(size(x,3),3);
    for i = 1:size(x,3)
        cx = convhull(x(:,:,i),'Simplify',true);
        [centrx,centry] = centroid(polyshape(x(cx(1:end-1),:,i)));
        r = max(sqrt( sum((x(cx(1:end-1),:)-[centrx,centry]).^2,2)));
        circs(:,i) = [ centrx, centry, r];
    end
end

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

function c = env_collision_avoidance(X,CONSTANTS)

    c = [];
    if isempty(CONSTANTS.obstacles) 
        return
    end
    c = zeros(2,size(CONSTANTS.obstacles,3));
    for i = 1:size(CONSTANTS.obstacles,3)
        obstacle = CONSTANTS.obstacles(:,:,i);
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
%             [~,simplex] = gjk2(p1,p2);
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

function c = vehicle_collision_avoidance(X1,X2,CONSTANTS)
    c = CONSTANTS.min_dist_intervehicles - MinDistBernstein2Polygon((X2-X1).',[0;0]);
end

function c = vehicle_collision_avoidance_bad(X1,X2,CONSTANTS)
    c = CONSTANTS.min_dist_intervehicles - min(sqrt(sum(BernsteinPow(X2-X1,2),2)));
end

function c = vehicle_collision_avoidance_isaac(X1,X2,CONSTANTS)
    c = CONSTANTS.min_dist_intervehicles -  min(sqrt(sum((CONSTANTS.BigElevMat*(X1-X2)).^2,2)));
end

function [c,ceq] = nonlcon(X,CONSTANTS)

    X = matrify(X,CONSTANTS);

    c = [];
    ceq = [];

    % dynamics
    ceq_dyn = cell(CONSTANTS.Nv,1);
    c_dyn = cell(CONSTANTS.Nv,1);
    for i = 1:CONSTANTS.Nv
        [c_i,ceq_i] = dynamics5vars(X(:,:,i),CONSTANTS);
        ceq_dyn{i} = ceq_i;
        c_dyn{i} = c_i;
    end
    ceq = [ceq; cell2mat(ceq_dyn)];
    c = [c; cell2mat(c_dyn)];
    
    % obstacles
    if ~isempty(CONSTANTS.obstacles) && ~CONSTANTS.uselogbar
        c_env_avoidance = cell(CONSTANTS.Nv,1);
        for i = 1:CONSTANTS.Nv
            c_env_avoidance{i} = env_collision_avoidance(X(:,:,i),CONSTANTS);
        end    
        c = [c; cell2mat(c_env_avoidance)];
    end
    if ~isempty(CONSTANTS.obstacles_circles) && ~CONSTANTS.uselogbar
        c_env_avoidance = cell(CONSTANTS.Nv,1);
        for i = 1:CONSTANTS.Nv
            c_env_avoidance{i} = env_collision_avoidance_circles(X(:,:,i),CONSTANTS);
        end
        c = [c; cell2mat(c_env_avoidance)];
    end

    % inter vehicles
    if CONSTANTS.Nv>1 && ~CONSTANTS.uselogbar
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

function J = costfun_single(X,CONSTANTS)
    v = X(:,4);
    w = X(:,5);
    a = CONSTANTS.DiffMat*v;
    J = sum(a.^2)+2*sum(w.^2);
    
end

function J = costfun(X,CONSTANTS)

    X = matrify(X,CONSTANTS);
    J = 0;

    % dynamics
    J_dyn = zeros(CONSTANTS.Nv,1);
    for i = 1:CONSTANTS.Nv
        J_dyn(i) = costfun_single(X(:,:,i),CONSTANTS);
    end    
    J = J + sum(J_dyn);

    % obstacles
    if ~isempty(CONSTANTS.obstacles) && CONSTANTS.uselogbar
        J_obs = cell(CONSTANTS.Nv,1);
        for i = 1:CONSTANTS.Nv
            J_obs{i} = arrayfun(@(z)logbarrierfunc(0.1,z), -env_collision_avoidance(X(:,:,i),CONSTANTS));
        end
        J = J+sum(cell2mat(J_obs));
    end
    if ~isempty(CONSTANTS.obstacles_circles) && CONSTANTS.uselogbar
        J_obs = cell(CONSTANTS.Nv,1);
        for i = 1:CONSTANTS.Nv
            J_obs{i} = arrayfun(@(z)logbarrierfunc(0.1,z), -env_collision_avoidance_circles(X(:,:,i),CONSTANTS));
        end
        J = J+sum(cell2mat(J_obs));

    end
    
    % inter vehicles
    if CONSTANTS.Nv>1  && true
        J_barrier = cell(nchoosek(CONSTANTS.Nv,2),1);
        it = 1;
        for i = 1:CONSTANTS.Nv
            for j = i+1:CONSTANTS.Nv
                J_barrier{it} = arrayfun(@(z)logbarrierfunc(0.5,z),-vehicle_collision_avoidance_isaac(X(:,1:2,i),X(:,1:2,j),CONSTANTS));          
                it = it+1;
            end
        end
        J = J+sum(cell2mat(J_barrier));
    end

end

function xinit = init_guess(CONSTANTS)

    N = CONSTANTS.N; 

    xinit = zeros(CONSTANTS.N-1,CONSTANTS.numvars,CONSTANTS.Nv);
    for i = 1:CONSTANTS.Nv
        x = linspace(CONSTANTS.xi(i,1),CONSTANTS.xf(i,1),N-1).';
        y = linspace(CONSTANTS.xi(i,2),CONSTANTS.xf(i,2),N-1).';
        %psi = -ones(N-1,1);
        psi = ones(N-1,1).*atan2(y(end)-y(1),x(end)-x(1));
        v = ones(N-1,1);
        w = zeros(N-1,1);
        xinit(:,:,i) = [x,y,psi,v,w];
    end

end

function xinit = bad_init_guess(CONSTANTS) %#ok<*DEFNU>

    N = CONSTANTS.N; 

    x = rand(N-1,1);
    y = rand(N-1,1);
    psi = rand(N-1,1);
    v = rand(N-1,1);
    w = rand(N-1,1);

    xinit = [x;y;psi;v;w];
    
end

function X = matrify(X,CONSTANTS)
    %X = [CONSTANTS.xi; reshape(X,[],CONSTANTS.numvars); CONSTANTS.xf];
    X = cat(1,...
        reshape(CONSTANTS.xi.',[1,CONSTANTS.numvars,CONSTANTS.Nv]),...
        reshape(X,[CONSTANTS.N-1,CONSTANTS.numvars,CONSTANTS.Nv]),...
        reshape(CONSTANTS.xf.',[1,CONSTANTS.numvars,CONSTANTS.Nv]));
end