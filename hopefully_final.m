clear; close all;

addpath('Bernstein');

%% Settings

N = 20; % order
T = 10; % time interval

%% Boundary Conditions
psi_0 = 0;
psi_f= 0 ;

init_conds = [3 5 cos(psi_0) sin(psi_0) 1]; % x y r1 r2 v
final_conds = [0 0 cos(psi_f) sin(psi_f) 0]; % inf <=> don't care


%% Linear Constraints

[Aeq,beq,lb,ub] = LinConstr(N,T,init_conds,final_conds);

%% Initial Guess

X = InitialGuess(N, T, init_conds, final_conds);

[~, ceq] = basicnonlcon(X,N,T);


%% DO the bloody thing

A = [];
b = [];

nonlcon=@(x) basicnonlcon(x,N,T);

options = optimoptions(@fmincon,'Algorithm','sqp',...
                       'MaxFunctionEvaluations',40000,...
                       'ConstraintTolerance',1e-6,...
                       'StepTolerance',1e-6,...
                       'Display','iter');
%                        'OptimalityTolerance',1e-8,...
% %                       'FunctionTolerance',1e-8,..
%                       'ConstraintTolerance',1e-3);

%options = optimoptions(@fmincon,'Algorithm','sqp');

tic;
    [xOut,Jout,exitflag,output] = fmincon(@(x)costFunc(x,N,T),X,A,b,Aeq,beq,lb,ub,nonlcon,options);
toc

%% Test the functions

axis equal
BernsteinPlot(xOut(:,1:2),T,'PlotControlPoints',false);

%% Cost function
function J = costFunc(X,N,T)

    [px, py, r1, r2, v, dx, dy, dr1, dr2, dv, ~, ~] = getControlPoints(X,N,T);
    rho = 100;

    deriv_mul = @(a,da,b,db) BernsteinSum(BernsteinMul(da,b),BernsteinMul(a,db));
    b_sq_int = @(p) BernsteinIntegr(BernsteinPow(p,2),T);

    % J1 = \int_0^T x^2 + dx^2 + rho ddx^2 dt
    J1 = b_sq_int(px) + b_sq_int(dx) + rho * b_sq_int(deriv_mul(v,dv,r1,dr1));
    J2 = b_sq_int(py) + b_sq_int(dy) + rho * b_sq_int(deriv_mul(v,dv,r2,dr2));
    J = norm([J1 J2]);

end


%% Linear contraints
% The linear contraints set the boundary conditions
function [Aeq, beq,lb,ub] = LinConstr(N,T,x0,xf)
    Aeq_boundary = zeros(5,(N+1)*7);

    Aeq_boundary(1,1)         = 1; % x0
    Aeq_boundary(2,1*(N+1)+1) = 1; % y0
    Aeq_boundary(3,2*(N+1)+1) = 1; % r10
    Aeq_boundary(4,3*(N+1)+1) = 1; % r20
    Aeq_boundary(5,4*(N+1)+1) = 1; % v0
    
    beq_boundary = x0';
    
    for i = 1:length(xf)
        if xf(i) ~= Inf
            Aeq_boundary(end+1, i*(N+1)) = 1;
            beq_boundary(end+1) = xf(i);
        end
    end

    Aeq_dyn = [ zeros(N+1,(N+1)*4), BernsteinDerivElevMat(N,T), zeros(N+1), -eye(N+1) ]; % dv = u2
    beq_dyn = zeros(N+1,1);

    Aeq = [Aeq_boundary; Aeq_dyn];
    beq = [beq_boundary; beq_dyn];
    
    % variable bounds
    
    % states
    x_ub = ones(5*(N+1),1)*Inf;
    x_lb = -x_ub;

    % control 
    u_ub = ones(2*(N+1),1)*Inf;
    u_lb = -u_ub; 

    lb = [x_lb(:); u_lb(:)];
    ub = [x_ub(:); u_ub(:)];

end

%% Nonlinear contraints
% The nonlinear contraints set the dynamic conditions
% this function is named "basic" because it doesn't take obstacles into
% account
function [c, ceq] = basicnonlcon(X,N,T)

    [~, ~, r1, r2, v, dx, dy, dr1, dr2, ~, u1, ~] = getTrajectories(X,N,T);

%     Mag = sqrt(r1.^2 + r2.^2);
%     r1 = (r1 + eps)./(Mag+eps); r2 = (r2 + eps)./(Mag + eps);
%     ceq = [
%         dx - v.*r1
%         dy - v.*r2
%         dr1 + r2.*u1
%         dr2 - r1.*u1
%     ]; 
    
    ceq = [
        dx - v.*r1
        dy - v.*r2
        dr1 + r2.*u1
        dr2 - r1.*u1
        r1.^2 + r2.^2 - 1
        %v.^2 - dx.^2 - dy.^2
    ];

    c = [];

end

%% Initial Guess
function [X] = InitialGuess(N, T, init_conds, final_conds)
    X = zeros(N+1,7);
    
    final_conds(final_conds==Inf) = 0;
    
    vi = init_conds(5)  .* init_conds(3:4);
    vf = final_conds(5) .* final_conds(3:4);
    
    X(1,1:5) = init_conds;
    X(end,1:5) = final_conds;
    
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
function [px, py, r1, r2, v, dx, dy, dr1, dr2, dv, u1, u2] = getControlPoints(X,N,T)
    
    persistent old_N old_T d_mat
    if  isempty(old_N) || isempty(old_T) || old_N ~= N || old_T ~= T
        old_N = N;
        old_T = T;
        d_mat = BernsteinDerivElevMat(N,T);
    end
    
    X = reshape(X(:), [N+1,7]); % this step could be redundant

    % get control points
    px = X(:,1);
    py = X(:,2);
    r1 = X(:,3);
    r2 = X(:,4);
    v = X(:,5);
    u1 = X(:,6);
    u2 = X(:,7);

    % compute  derivatives 
    dX  = d_mat*X;
    dx  = dX(:,1);
    dy  = dX(:,2);
    dr1 = dX(:,3);
    dr2 = dX(:,4);
    dv  = dX(:,5);

end

%% Get Control Points of Trajectories
function [px, py, r1, r2, v, dx, dy, dr1, dr2, dv, u1, u2] = getTrajectories(X,N,T)
    
    persistent old_N old_T d_mat eval_mat
    if  isempty(old_N) || isempty(old_T) || old_N ~= N || old_T ~= T
        old_N = N;
        old_T = T;
        d_mat = BernsteinDerivElevMat(N,T);
        eval_mat = BernsteinCtrlPntsEval(N);
    end
    
    X = reshape(X(:), [N+1,7]); % this step could be redundant

    dX = d_mat*X;
    X  = eval_mat*X;
    dX = eval_mat*dX;
    
    
    % get points
    px = X(:,1);
    py = X(:,2);
    r1 = X(:,3);
    r2 = X(:,4);
    v = X(:,5);
    u1 = X(:,6);
    u2 = X(:,7);

    % compute  derivatives 

    dx  = dX(:,1);
    dy  = dX(:,2);
    dr1 = dX(:,3);
    dr2 = dX(:,4);
    dv  = dX(:,5);

end