clear all; %#ok<CLALL>

addpath('..\Bernstein');
addpath('..\BeBOT_lib');
addpath('..\TrajecOptimLib');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

modelparameters = struct();
modelparameters.beta = 1;
constants.modelparameters = modelparameters;

%%%%%%%%%%%%%%%% Example 0 %%%%%%%%%%%%%%%%
constants.T = 100;
%                x  y  psi  u  v w
constants.xi = [ 0  0  pi/4 .1 0 0 ]; % initial conds
constants.xf = [ 10 10 pi/4 .1 0 0 ]; % final conds
constants.N = 20;
constants.obstacles_circles = [ 5, 5, 1];
constants.statebounds = [
    -Inf, -Inf, -2*pi, 0, -Inf, -1; % inferior bounds of each state variable
    Inf, Inf, 2*pi, 5, Inf, 1; % superior bounds of each state variable
];

constants.numinputs = 2;

%%%%%%%%%%%%%%%% Functions %%%%%%%%%%%%%%%%
constants.costfun_single = @costfun_single;
constants.dynamics = @dynamicsHovercraft;
% constants.init_guess = @init_guess;
constants.recoverxy = @recoverplot;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[xOut,JOut] = run_problem(constants); % runs once for a desired N
% [xOut,JOut] = run_problem_progressive_n(constants); % runs with increasing N

%%

disp(['The final cost is ', num2str(JOut)])
%constraint_evaluation(xOut,constants);
plot_xy(xOut,constants);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t,xy, tenpoints] = recoverplot(X,constants)

    tau_u = @(t) BernsteinEval(X(:,7),constants.T,t);
    tau_w = @(t) BernsteinEval(X(:,8),constants.T,t);

    [t,xy] = ode45(@(t,xy)odefunc(t,xy,tau_u, tau_w, constants.modelparameters), [0 constants.T], X(1,1:6));

    tenpoints =  BernsteinEval(X,constants.T,linspace(0,constants.T,10));

    function dydt = odefunc(t,y,tau_u,tau_w,modelparameters)
        
        yaw = y(3);
        u = y(4);
        v = y(5);
        w = y(6);
        
        dydt = zeros(6,1);
        dydt(1) = u*cos(yaw) - v*sin(yaw); %dx
        dydt(2) = u*sin(yaw) + v*cos(yaw); % dy
        dydt(3) = w; % dyaw
        dydt(4) = v*w + tau_u(t); %du
        dydt(5) = - u*w - modelparameters.beta*v; % dv
        dydt(6) = tau_w(t); % dr  i
        
    end

end


function [c,ceq] = dynamicsHovercraft(X,constants)

    DiffMat = constants.DiffMat;

    xp = X(:,1);
    yp = X(:,2);
    psi = X(:,3);
    u = X(:,4);
    v = X(:,5);
    w = X(:,6);
    tau_u = X(:, 7);
    tau_w = X(:, 8);

    ceq = [
        DiffMat*xp - v.*cos(psi);
        DiffMat*yp - v.*sin(psi);
        DiffMat*psi - w;
        DiffMat*u - v.*w - tau_u;
        DiffMat*v + u.*w + constants.modelparameters.beta*v;
        DiffMat*w - tau_w;
    ];

    c = [];

end


function J = costfun_single(X,constants)
    % the matrix X's columns have, control points for x,y,psi,v,w variables for 1 vehicle,
    % respectfully

    v = X(:,4);
    w = X(:,6);
    dv = constants.DiffMat*v;
    dw = constants.DiffMat*w;
    J = constants.T/(constants.N+1)*sum(2*dw.^2+dv.^2);

end


function xinit = init_guess(constants)

    N = constants.N; 

    xinit = zeros((constants.N-1)*constants.numvars,constants.Nv);
    for i = 1:constants.Nv
        x = linspace(constants.xi(i,1),constants.xf(i,1),N-1).';
        y = linspace(constants.xi(i,2),constants.xf(i,2),N-1).';
        psi = ones(N-1,1).*atan2(y(end)-y(1),x(end)-x(1));
        u = ones(N-1,1);
        v = ones(N-1,1);
        w = zeros(N-1,1);
        xinit(:,i) = [x;y;psi;u;v;w];
    end
    xinit = [xinit ; zeros((constants.N+1)*constants.numinputs,constants.Nv)];

end
