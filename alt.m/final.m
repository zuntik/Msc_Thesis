clear; close all;

addpath('Bernstein');

%% Settings

N = 40; % order
T = 10; % time interval

%% Boundary Conditions
psi_0 = 0;
psi_f= 0 ;

init_conds = [3 5 cos(psi_0) sin(psi_0) 1 0]; % x y r1 r2 v w
final_conds = [0 0 cos(psi_f) sin(psi_f) 0 0]; % inf <=> don't care


%% Linear Constraints

[Aeq,beq,lb,ub] = LinConstr(N,T,init_conds,final_conds);

%% Initial Guess

xIn = InitialGuess(N, T, init_conds, final_conds);

[~, ceq] = basicnonlcon(xIn,N,T);


%% DO the bloody thing

A = [];
b = [];

nonlcon=@(x) basicnonlcon(x,N,T);

% options = optimoptions(@fmincon,'Algorithm','sqp',...
%                        'MaxFunctionEvaluations',40000,...
%                        'ConstraintTolerance',1e-6,...
%                        'StepTolerance',1e-6,...
%                        'Display','iter');
% %                        'OptimalityTolerance',1e-8,...
% % %                       'FunctionTolerance',1e-8,..
% %                       'ConstraintTolerance',1e-3);

options = optimoptions(@fmincon,'Algorithm','sqp');

tic;
    [xOut,Jout,exitflag,output] = fmincon(@(x)costFunc(x,N,T),xIn,A,b,Aeq,beq,lb,ub,nonlcon,options);
toc

%% Test the functions

axis equal
BernsteinPlot(xOut(:,1:2),T,'PlotControlPoints',false);

%% Cost function
function J = costFunc(X,N,T)

    [~, ~, ~, ~, ~, ~, ~, ~, v, w]  = getControlPoints(X,N,T);

    b_sq_int = @(p) BernsteinIntegr(BernsteinPow(p,2),T);

    cw = 1;
    cv = 1;
    
    J = cw * b_sq_int(w) + cv* b_sq_int(v);

end


%% Linear contraints
% The linear contraints set the boundary conditions
function [Aeq, beq,lb,ub] = LinConstr(N,T,xi,xf)
    
    Aeq = zeros(6,(N+1)*6);

    for i = 1:6
        Aeq(i,1+(i-1)*(N+1)) = 1;
    end

    beq = xi';

    for i = 1:length(xf)
        if xf(i) ~= Inf
            Aeq(end+1, i*(N+1)) = 1;
            beq(end+1) = xf(i);
        end
    end

    % variable bounds
    
    % states
    x_ub = ones(4*(N+1),1)*Inf;
    x_lb = -x_ub;

    % control 
    
    vmax = norm(xf(1:2)-xi(1:2))/T*1.5;
    u_ub = ones((N+1),2).*[vmax 10*180/pi ];
    u_lb = ones((N+1),2).*[0 -10*180/pi]; 

    lb = [x_lb(:); u_lb(:)];
    ub = [x_ub(:); u_ub(:)];

end

%% Nonlinear contraints
% The nonlinear contraints set the dynamic conditions
% this function is named "basic" because it doesn't take obstacles into
% account
function [c, ceq] = basicnonlcon(X,N,T)

    [~, ~, c, s, dx, dy, dc, ds, w, v] = getTrajectories(X,N,T);

    ceq = [
        dx - v.*c
        dy - v.*s
        dc + s.*w
        ds - c.*w
        c.^2 + s.^2 - 1
        %v.^2 - dx.^2 - dy.^2
    ];

    c = [];

end

%% Initial Guess
function [X] = InitialGuess(N, T, init_conds, final_conds)
    X = zeros(N+1,6);
    
    final_conds(final_conds==Inf) = 0;
    
    vi = init_conds(5)  .* init_conds(3:4);
    vf = final_conds(5) .* final_conds(3:4);
    
    X(1,1:6) = init_conds;
    X(end,1:6) = final_conds;
    
    min_dist = norm(final_conds(1:2)-init_conds(1:2));
    
    % deal with the second and penultimate control points
    X(2,1:2) = (T/N).*vi + init_conds(1:2);
    X(end-1,1:2) = -(T/N).*vf + final_conds(1:2);
    
    % deal with third and antepenultimate
    X(3,1:2) =     X(2,1:2)     + (min_dist/N).* init_conds(3:4);
    X(end-2,1:2) = X(end-1,1:2) - (min_dist/N).*final_conds(3:4);
    
    % place remaning control points in a straight line
    X(3:end-2,1) = linspace(X(3,1),X(end-2,1),N-3);
    X(3:end-2,2) = linspace(X(3,2),X(end-2,2),N-3);
    
    vec = X(end-2,1:2) - X(3,1:2);
    
    psi_mid = atan2(vec(2),vec(1));

    for it = 3:(size(X,1)-3)
        X(it,3:4) = [ cos(psi_mid) sin(psi_mid) ];
    end
    X(2,3:4) = init_conds(3:4);
    vec_fin = X(end,1:2)-X(end-2,1:2);
    for it = (-2:0)+size(X,1)
        psi_fin = atan2(vec_fin(2),vec_fin(1));
        X(it,3:4) = [ cos(psi_fin) sin(psi_fin) ];
    end
    
end


%% Get Control Points of Trajectories
function [px, py, c, s, dx, dy, dc, ds, v, w] = getControlPoints(X,N,T)
    
    persistent old_N old_T d_mat
    if  isempty(old_N) || isempty(old_T) || old_N ~= N || old_T ~= T
        old_N = N;
        old_T = T;
        d_mat = BernsteinDerivElevMat(N,T);
    end
    
    X = reshape(X(:), [N+1,6]); % this step could be redundant

    % get control points
    px = X(:,1);
    py = X(:,2);
    c = X(:,3);
    s = X(:,4);
    v = X(:,5);
    w = X(:,6);

    % compute  derivatives 
    dX  = d_mat*X;
    dx  = dX(:,1);
    dy  = dX(:,2);
    dc = dX(:,3);
    ds = dX(:,4);

end

%% Get Trajectories
function [px, py, c, s, dx, dy, dc, ds, w, v] = getTrajectories(X,N,T)
    
    persistent old_N old_T d_mat eval_mat
    if  isempty(old_N) || isempty(old_T) || old_N ~= N || old_T ~= T
        old_N = N;
        old_T = T;
        d_mat = BernsteinDerivElevMat(N,T);
        eval_mat = BernsteinCtrlPntsEval(N);
    end
    
    X = reshape(X(:), [N+1,6]); % this step could be redundant

    dX = d_mat*X;
    X  = eval_mat*X;
    dX = eval_mat*dX;
    
    % get control points
    px = X(:,1);
    py = X(:,2);
    c = X(:,3);
    s = X(:,4);
    v = X(:,5);
    w = X(:,6);

    % compute  derivatives 
    dx  = dX(:,1);
    dy  = dX(:,2);
    dc = dX(:,3);
    ds = dX(:,4);

end

