clear all; close all; %#ok<CLALL>

addpath('..\Bernstein');
addpath('..\BeBOT_lib');
addpath('..\epa');

% Settings

CONSTANTS.T = 10; % time interval
% CONSTANTS.T = 50; % time interval
CONSTANTS.xi = [
    5 3 0 1 0
    0 1 0 1 0
    ]; % x y psi v w
CONSTANTS.xf = [
    0 0 0 1 0
    4 1 0 1 0 
    ]; % x y psi v w
% CONSTANTS.xi  = [-10 40 0 1.1 0]; % x y psi v w
% CONSTANTS.xf = [0 0 0 1.2 0]; % x y psi v w
CONSTANTS.N = 40;
CONSTANTS.min_dist_intervehicles = 0.1;
CONSTANTS.DiffMat = BernsteinDerivElevMat(CONSTANTS.N,CONSTANTS.T);
CONSTANTS.EvalMat = BernsteinCtrlPntsEval(CONSTANTS.N);
CONSTANTS.BigElevMat = BernsteinDegrElevMat(CONSTANTS.N,CONSTANTS.N*10);
CONSTANTS.numvars = size(CONSTANTS.xi,2);
CONSTANTS.Nv = size(CONSTANTS.xi,1);%number of vehicles
% CONSTANTS.obstacles = [  1 1 ; 1 2 ; 2 2 ; 2 1 ] + [100 100];
CONSTANTS.obstacles = [];

xIn = init_guess(CONSTANTS);

options = optimoptions(@fmincon,'Display','Iter','Algorithm','sqp','MaxFunctionEvaluations',300000,'StepTolerance',eps,'MaxIterations',Inf);

xOut_uc = fmincon(@(x)costfun(x,CONSTANTS),xIn,[],[],[],[],[],[],@(x)nonlcon(x,CONSTANTS),options);
%CONSTANTS.obstacles = [  1 1 ; 1 2 ; 2 2 ; 2 1 ];%+[10 10];
if ~isempty(CONSTANTS.obstacles)
    xOut = fmincon(@(x)costfun(x,CONSTANTS),xOut_uc,[],[],[],[],[],[],@(x)nonlcon(x,CONSTANTS),options);
else
    xOut = xOut_uc;
end

%%

[c,ceq] = nonlcon(xOut,CONSTANTS);
disp((c<0).')
disp(norm(ceq))
disp(ceq.')

%%

xOut = [CONSTANTS.xi; reshape(xOut,[],CONSTANTS.numvars); CONSTANTS.xf];
%%
figure
BernsteinPlot(xOut(:,1:2),CONSTANTS.T);
[~,xy] = recoverplot(xOut,CONSTANTS.T);
plot(xy(:,1),xy(:,2));

if ~isempty(CONSTANTS.obstacles)
    for i = 1:size(CONSTANTS.obstacles,3)
        plot(polyshape(CONSTANTS.obstacles(:,:,i)))
    end
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
    c = [-v+0.2;psi-pi;-psi-pi];
    %c = [ -v; psi-pi-pi/3; -psi-pi-pi/3];
    
end

function c = env_collision_avoindance(X,CONSTANTS)
    
    if isempty(CONSTANTS.obstacles)
        c = [];
        return
    end
    
    N = CONSTANTS.N;
    
    cpts = X(:,1:2);
    curvepoints = CONSTANTS.EvalMat*X(:,1:2);

    c = -1.*ones(N,size(CONSTANTS.obstacles,3));
    for it = 1:size(CONSTANTS.obstacles,3)
        obstacle = CONSTANTS.obstacles(:,:,it);
        pobs = ConvexPolygon(obstacle);
        intersecting_indexes = recursive_collision_test(cpts,curvepoints,1,N+1,pobs);
        if ~any(intersecting_indexes)
            continue
        end
        for i = 1:N
            if intersecting_indexes(i)&&intersecting_indexes(i+1)
                pcurv = ConvexPolygon([ cpts(i:i+1,:); curvepoints(i:i+1,:)]);
                [~,simplex] = gjk2(pobs,pcurv);
                support=@(d) pobs.support(d)-pcurv.support(-d);
                vec_up = directed_epa(support,simplex,[0;1]);
                c(i,it) = vec_up(2);
            end
        end
    end
    c = c(:);
    
end

function c = env_collision_avoindance_good(X,CONSTANTS)

    c = [-1;-1 ];
    
    if isempty(CONSTANTS.obstacles)
        return
    end
    
    N = CONSTANTS.N;
    
    cpts = X(:,1:2);
    curvepoints = CONSTANTS.EvalMat*X(:,1:2);

    for it = 1:size(CONSTANTS.obstacles,3)
        obstacle = CONSTANTS.obstacles(:,:,it);
        pobs = ConvexPolygon(obstacle);    
        intersecting_indexes = recursive_collision_test(cpts,curvepoints,1,N,pobs);
        if ~any(intersecting_indexes)
            continue
        end
        vertices = [ cpts(intersecting_indexes,:); curvepoints(intersecting_indexes,:)];
        pcurv = ConvexPolygon(vertices);
        [~,simplex] = gjk2(pobs,pcurv);
        support=@(d) pobs.support(d)-pcurv.support(-d);
        vec_up = directed_epa(support, simplex, [0;1]);
        %vec_side = directed_epa(support, simplex, [1;0]);
    end
    
    if ~any(intersecting_indexes)
        %c = [-1;-1];
        c  = -1;
    else
        %c = [vec_up(2); vec_side(1)];
        c = vec_up(2);
    end
    
end

function indexes = recursive_collision_test(ctrl_pts,curv_pts,l,r,p_obs)

    p_curv = ConvexPolygon([ctrl_pts(l:r,:);curv_pts(l:r,:)]);
    indexes = false(r-l+1,1);
    if ~gjk2(p_obs,p_curv)
        return
    else
        if r-l == 1
            indexes = true(2,1);
            return
        else
            halfsize = floor((r-l)/2)+1;
            midp = l+halfsize-1;% middle point
            indexes(1:halfsize) = indexes(1:halfsize)|recursive_collision_test(ctrl_pts,curv_pts,l,midp,p_obs);
            indexes(halfsize:end) = indexes(halfsize:end)|recursive_collision_test(ctrl_pts,curv_pts,midp,r,p_obs);
        end
    end

end

function c = env_collision_avoindance_fair(X,CONSTANTS)

    c = [-1;-1 ];
    
    if isempty(CONSTANTS.obstacles)
        return
    end
    
    N = CONSTANTS.N;
    
    cpts = X(:,1:2);
    curvepoints = CONSTANTS.EvalMat*X(:,1:2);

    for it = 1:size(CONSTANTS.obstacles,3)
        obstacle = CONSTANTS.obstacles(:,:,it);
        p1 = ConvexPolygon(obstacle);    
        allpoligons = cell(N,1);
        intersecting_indexes = false(N,1);
        for i = 1:N
            if i == 1
                vertices = [cpts(i:i+1,:);curvepoints(i+1,:)];
            elseif i == N
                vertices = [cpts(i:i+1,:);curvepoints(i,:)];
            else
                vertices = [cpts(i:i+1,:);curvepoints(i:i+1,:)];
            end
            vertices = uniquetol(vertices,'ByRows',true);
            allpoligons{i} = vertices;
            p2 = ConvexPolygon(vertices);
            intersecting_indexes(i) = gjk2(p1,p2);
        end
        if ~any(intersecting_indexes)
            continue
        end
        vertices = cell2mat(allpoligons(intersecting_indexes));
        p2 = ConvexPolygon(uniquetol(vertices,'ByRows',true));
        [~,simplex] = gjk2(p1,p2);
        support=@(d) p1.support(d)-p2.support(-d);
        vec_up = directed_epa(support, simplex, [0;1]);
        vec_side = directed_epa(support, simplex, [1;0]);
    end
    
    if ~any(intersecting_indexes)
        c = [-1;-1];
    else
        c = [vec_up(2); vec_side(1)];
    end
    
end

function c = env_collision_avoindance_bad(X,CONSTANTS)

    c = [];
    if ~isempty(CONSTANTS.obstacles) 
        for i = 1:size(CONSTANTS.obstacles,3)
            obstacle = CONSTANTS.obstacles(:,:,i);
            [dist, t, pt] = MinDistBernstein2Polygon_extended(X(:,1:2).', obstacle.');
            if dist>0
                c = [c; -1; -1]; %#ok<AGROW>
            else
                if length(t) == 2 && false
                    cptsleftcut = deCasteljau(X(:,1:2).',t(1));
                    cptsrightcut = deCasteljau(cptsleftcut(:,N+1:end),t(2));
                    cpts_inside = cptsrightcut(:,1:N+1);
                else
                    cpts_inside = pt.';
                end
                p2 = ConvexPolygon(cpts_inside.');
                p1 = ConvexPolygon(obstacle);
                [~,simplex] = gjk2(p1,p2);
                support=@(d) p1.support(d)-p2.support(-d);
                vec_up = directed_epa(support, simplex, [0;1]);
                vec_side = directed_epa(support, simplex, [1;0]);
                c = [c; vec_up(2); vec_side(1)]; %#ok<AGROW>
            end
        end
    end
    
end

function c = env_collision_avoindance_more_or_less(X,CONSTANTS)

    c = [-1;-1 ];
    
    if isempty(CONSTANTS.obstacles)
        return
    end
    
    N = CONSTANTS.N;
    
    cpts = X(:,1:2);
    curvepoints = CONSTANTS.EvalMat*X(:,1:2);

    for it = 1:size(CONSTANTS.obstacles,3)
        obstacle = CONSTANTS.obstacles(:,:,it);
        p1 = ConvexPolygon(obstacle);    
        upmin = -1;
        sidemin = -1;

        allpoligons = cell(N);
        for i = 1:N
            if i == 1
                vertices = [cpts(i:i+1,:);curvepoints(i+1,:)];
            elseif i == N
                vertices = [cpts(i:i+1,:);curvepoints(i,:)];
            else
                vertices = [cpts(i:i+1,:);curvepoints(i:i+1,:)];
            end
            vertices = uniquetol(vertices,'ByRows',true);
            allpoligons{i} = sortvertices(vertices);
            p2 = ConvexPolygon(vertices);
            [intersection,simplex] = gjk2(p1,p2);
            if intersection
                support=@(d) p1.support(d)-p2.support(-d);
                vec_up = directed_epa(support, simplex, [0;1]);
                vec_side = directed_epa(support, simplex, [1;0]);
                if vec_up(2) > upmin 
                    upmin = vec_up(2);
                end
                if vec_side(1) > sidemin
                    sidemin = vec_side(1);
                end
            end
        end
    end
    
    c = [upmin; sidemin];
    
end

function c = vehicle_avoidance(X1,X2)
    X = X2-X1;
    X = CONSTANTS.BigElevMat * X;
    mag_squared = sum(X.^2,2);
    c = CONSTANTS.min_dist_intervehicles-mag_squared;
end

function [c,ceq] = nonlcon(X,CONSTANTS)

    X = cat(1,...
        reshape(CONSTANTS.xi.',[1,CONSTANTS.numvars,CONSTANTS.Nv]),...
        reshape(X,[CONSTANTS.N+1,CONSTANTS.numvars,CONSTANTS.Nv]),...
        reshape(CONSTANTS.xi.',[1,CONSTANTS.numvars,CONSTANTS.Nv]));
    
    ceq = cell(CONSTANTS.Nv,1);
    c = cell(CONSTANTS.Nv,1);
    c_avoidance = cell(CONSTANTS.Nv,1);
    
    for i = 1:CONSTANTS.Nv
        [c_i,ceq_i] = dynamics5vars(X(:,:,i),CONSTANTS);
        c_i_avoidance = env_collision_avoindance(X(:,:,i),CONSTANTS);
        ceq{i} = ceq_i;
        c{i} = c_i;
        c_avoidance{i} = c_i_avoidance;
    end
    ceq=cell2mat(ceq);
    c = cell2mat(c);
    c_avoidance = cell2mat(c_avoidance);
    
    % only perform inter vehicular collision detection if Nv>1
    if CONSTANTS.Nv > 1
        c_vehicle_avoidance = cell(nchoosek(CONSTANTS.Nv,2));
        maycollide = 1;
    else
        c_vehicle_avoidance = [];
        maycollide = 0;
    end
    
    it = 1;
    for i = 1:CONSTANTS.Nv*maycollide
        for j = i+1:Nv
            c_vehicle_avoidance{it} = vehicle_avoidance(X(:,1:2,i),X(:,1:2,j));
            it = it+1;
        end
    end
    c_vehicle_avoidance = cell2mat(c_vehicle_avoidance);
        
    c = [c;c_avoidance;c_vehicle_avoidance];
    
end

function J = costfun(X,CONSTANTS)

    X = [CONSTANTS.xi; reshape(X,[],CONSTANTS.numvars); CONSTANTS.xf];
    v = X(:,4);
    w = X(:,5);
    a = CONSTANTS.DiffMat*v;
    J = sum(a.^2)+2*sum(w.^2);
    
end

function xinit = init_guess(CONSTANTS)

    N = CONSTANTS.N; 

    x = linspace(CONSTANTS.xi(1),CONSTANTS.xf(1),N-1).';
    y = linspace(CONSTANTS.xi(2),CONSTANTS.xf(2),N-1).';
    psi = -ones(N-1,1);
    v = ones(N-1,1);
    w = zeros(N-1,1);

    xinit = [x;y;psi;v;w];
    
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
