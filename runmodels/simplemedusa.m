clear all; %#ok<CLALL>

addpath('..\Bernstein');
addpath('..\BeBOT_lib');
addpath('..\TrajecOptimLib');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% medusa model parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MODELPARAMS = struct();

% %%%%%% Coefficients
% MODELPARAMS.mass = 17.0;
% MODELPARAMS.I_z = 4.14;
% MODELPARAMS.X_dot_u = -20;
% MODELPARAMS.Y_dot_v = -1.3175;
% MODELPARAMS.N_dot_r =  -8.69;
% MODELPARAMS.X_u = -0.2;
% MODELPARAMS.Y_v = -55.117;
% MODELPARAMS.N_r = -4.14;
% MODELPARAMS.X_uu = -25; 
% MODELPARAMS.Y_vv = -101.2776;
% MODELPARAMS.N_rr = -6.23;

%%%%%% Coefficients
MODELPARAMS.mass = 17.0;
MODELPARAMS.I_z = 1;
MODELPARAMS.X_dot_u = -20;
MODELPARAMS.Y_dot_v = -30;%-1.3175;
MODELPARAMS.N_dot_r = -8.69;% -0.5;
MODELPARAMS.X_u = -0.2;
MODELPARAMS.Y_v = -50;
MODELPARAMS.N_r = -4.14; %-0.1
MODELPARAMS.X_uu = -25;
MODELPARAMS.Y_vv = -0.01;%-101.2776;
MODELPARAMS.N_rr = -6.23; %-21

%%%%%% masses
MODELPARAMS.m_u = MODELPARAMS.mass - MODELPARAMS.X_dot_u;
MODELPARAMS.m_v = MODELPARAMS.mass - MODELPARAMS.Y_dot_v;
MODELPARAMS.m_r = MODELPARAMS.I_z - MODELPARAMS.N_dot_r;
MODELPARAMS.m_uv = MODELPARAMS.m_u - MODELPARAMS.m_v;

%%%%%% Constants
MODELPARAMS.fu  = 0;
MODELPARAMS.fv  = 0;
MODELPARAMS.fr  = 0;
MODELPARAMS.Vcx = 0;
MODELPARAMS.Vcy = 0;

constants.MODELPARAMS=MODELPARAMS;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% runtime parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% % 1 vehicle no constraints
% constants.T = 80; % time
% constants.xi = [0 0 0    .7 0 0]; % x y yaw u v r
% constants.xf = [40 2 pi/2 .7 0 0];
% constants.plotboatsize = 1;
% % constants.obstacles_circles = [8, 0, 2];
% constants.N = 15;

constants.T = 100;
constants.xi = [
    0 -5 0 .9 0 0
    0 5 0 .9 0 0
    0 0 0 .9 0 0
    ];
constants.xf = [
    40 5 0 .9 0 0
    40 -5 0 .9 0 0
    40 0 0 .9 0 0
    ];
constants.N = 10;

% constants.T = 60; % time
% constants.xi = [0 0 0    1 0 0]; % x y yaw u v r
% constants.xf = [30 30 pi/2 1 0 0]; % x y yaw u v r
% constants.xf = [30 5 pi/6 1 0 0];
% constants.xf = [ 0 6 pi 1 0 0 ];

% 3 vehicles 1 circle obstacle
% constants.N = 40;
% constants.T = 20;
% constants.xi = [
%     -10 4 0 1 0 0
% %     -10 -4 0 1 0 0
% %     -10 0 0 1 0 0
% ];
% constants.xf = [
%     10 -1 0 1 0 0
% %     10 1 0 1 0 0
% %     10 0 0 1 0 0
% ];
% constants.T = 15;
% constants.obstacles = [];
% %constants.obstacles_circles = [ 0 0 3]; % x y r
% constants.obstacles_circles = [];


% common parameters
constants.min_dist_int_veh = 3;
constants.numinputs = 2;
% constants.uselogbar = false; % use log barrier func
% constants.usesigma = true; % a variable for the usage of log barrier func

constants.statebounds = [
    % x   y    yaw   u   v     r
    -1000 -1000 -1000 0 -1000 -.74 % lower state bounds
    1000 1000 1000 1.1 1000 .74 % upper state bounds
    ];
constants.inputbounds = [
    % t_u t_r
    0 -.113; % lower input bounds
    25.9 .113 % upper input bounds
    ];


%%% bad constraints
% constants.plotboatsize = 0.5;
% constants.T = 15; % time
% constants.xi = [0 0 0    1 0 0]; % x y yaw u v r
% constants.xf = [5 5 pi/2 1 0 0]; % x y yaw u v r
% % constants.xf = [5 2 pi/6 1 0 0];
% % constants.xf = [ 0 6 pi 1 0 0 ];
% % constants.xf = [15 2 pi/2 1 0 0];
% constants.statebounds = [
%   %     %   y    yaw   u   v     r
%     -1000 -1000 -1000 0 -1000 -1000 % lower state bounds
%     1000 1000 1000 1000 1000 1000 % upper state bounds
%     ];
% constants.inputbounds = [
%     % t_u t_r
%     0 -5 % lower input bounds
%     1000 5 % upper input bounds
%     ];  


% functions
constants.costfun_single = @costfun_single;
constants.dynamics = @dynamicsmedusa;
constants.init_guess = @init_guess;
constants.recoverxy = @recoverplot;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xOut = run_problem_progressive_n(constants);
%xOut = run_problem(constants);

%%
constants = processconstants(constants);

%constraint_evaluation(xOut,constants);
plot_xy(xOut,constants);

points = BernsteinEval(xOut,constants.T,linspace(0,constants.T,10));
for v = 1:size(constants.xi, 1)
    for i = 1:10
        plotboat(points(i,1,v),points(i,2,v),points(i,3,v),constants.plotboatsize);
    end
end

if size(constants.xi,1) == 2
    figure
    X_diff = xOut(:,1:2,1)-xOut(:,1:2,2);
    X_diff_magsquared = sum(BernsteinPow(X_diff,2),2);
%     fplot(@(t)sqrt(BernsteinEval(X_diff_magsquared,constants.T,t)),[0, constants.T]);
    fplot(@(t)sqrt(sum(BernsteinEval(X_diff,constants.T,t).^2,2)), [0,constants.T]);
%     BernsteinPlot(Mag,constants.T);
end

if isfield(constants, 'obstacles_circles') && ~isempty(constants.obstacles_circles) &&...
        size(constants.obstacles_circles,3) == 1 && constants.Nv == 1
    figure, grid on
    BernsteinPlot(sum((xOut(:,1:2)-constants.obstacles_circles(1:2)).^2,2),constants.T);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t,xy] = recoverplot(X,constants)

    tau_u = @(t) BernsteinEval(X(:,7),constants.T,t);
    tau_r = @(t) BernsteinEval(X(:,8),constants.T,t);
    
    [t,xy] = ode45(@(t,xy)odefunc(t,xy,tau_u,tau_r,constants.MODELPARAMS), [0 constants.T], X(1,1:6));

    function dydt = odefunc(t,y,tau_u,tau_r,MODELPARAMS)

        % load the parameters
        %%%%%% Coefficients
%         mass = MODELPARAMS.mass;
%         I_z = MODELPARAMS.I_z;
%         X_dot_u = MODELPARAMS.X_dot_u;
%         Y_dot_v = MODELPARAMS.Y_dot_v;
%         N_dot_r = MODELPARAMS.N_dot_r;
        X_u = MODELPARAMS.X_u;
        Y_v = MODELPARAMS.Y_v;
        N_r = MODELPARAMS.N_r;
        X_uu = MODELPARAMS.X_uu;
        Y_vv = MODELPARAMS.Y_vv;
        N_rr = MODELPARAMS.N_rr;

        %%%%%% masses
        m_u = MODELPARAMS.m_u;
        m_v = MODELPARAMS.m_v;
        m_r = MODELPARAMS.m_r;
        m_uv = MODELPARAMS.m_uv;

        %%%%%% Constants
        fu = MODELPARAMS.fu;
        fv = MODELPARAMS.fv;
        fr = MODELPARAMS.fr;
        Vcx = MODELPARAMS.Vcx;
        Vcy = MODELPARAMS.Vcy;

        yaw = y(3);
        u = y(4);
        v = y(5);
        r = y(6);

        %%%%%%%%% drag
        d_u = -X_u - X_uu*abs(u);
        d_v = -Y_v - Y_vv*abs(v);
        d_r = -N_r - N_rr*abs(r);

        dydt = zeros(6,1);
        dydt(1) = u*cos(yaw) - v*sin(yaw) + Vcx; %dx
        dydt(2) = u*sin(yaw) + v*cos(yaw) + Vcy; % dy
        dydt(3) = r; % dyaw
        dydt(4) = 1/m_u*(tau_u(t) + m_v*v*r - d_u*u+fu); %du
        dydt(5) = 1/m_v*(-m_u*u*r - d_v*v+fv); % dv
        dydt(6) = 1/m_r*(tau_r(t) + m_uv*u*v - d_r*r+fr); % dr  
    end

end

function [t,xy] = recoverplot_simulink(X,constants)

    h=get_param('simplemedusa_model','modelworkspace');
    h.assignin('constants',constants);
    h.assignin('init_vals', X(1,1:6));
    h.assignin('tau_u', X(:,7));
    h.assignin('tau_r', X(:,8));
    out = sim('simplemedusa_model', constants.T);
    t = out.x.time;
    xy = [ 
        out.x.signals.values, ...
        out.y.signals.values, ...
        out.yaw.signals.values, ...
        out.u.signals.values, ...
        out.v.signals.values, ...
        out.r.signals.values
        ];
end

function [c,ceq] = dynamicsmedusa(X,constants)

    DiffMat = constants.DiffMat;
    
    % states
    x = X(:,1);
    y = X(:,2);
    yaw = X(:,3);
    u = X(:,4);
    v = X(:,5);
    r = X(:,6);
    % inputs
    tau_u = X(:,7);
    tau_r = X(:,8); 

    % load the parameters
    %%%%%% Coefficients
%     mass = constants.MODELPARAMS.mass;
%     I_z = constants.MODELPARAMS.I_z;
%     X_dot_u = constants.MODELPARAMS.X_dot_u;
%     Y_dot_v = constants.MODELPARAMS.Y_dot_v;
%     N_dot_r = constants.MODELPARAMS.N_dot_r;
    X_u = constants.MODELPARAMS.X_u;
    Y_v = constants.MODELPARAMS.Y_v;
    N_r = constants.MODELPARAMS.N_r;
    X_uu = constants.MODELPARAMS.X_uu;
    Y_vv = constants.MODELPARAMS.Y_vv;
    N_rr = constants.MODELPARAMS.N_rr;

    %%%%%% masses
    m_u = constants.MODELPARAMS.m_u;
    m_v = constants.MODELPARAMS.m_v;
    m_r = constants.MODELPARAMS.m_r;
    m_uv = constants.MODELPARAMS.m_uv;

    %%%%%% Constants
    fu = constants.MODELPARAMS.fu;
    fv = constants.MODELPARAMS.fv;
    fr = constants.MODELPARAMS.fr;
    Vcx = constants.MODELPARAMS.Vcx;
    Vcy = constants.MODELPARAMS.Vcy;

    %%%%%%%%% drag
    d_u = -X_u - X_uu*abs(u);
    d_v = -Y_v - Y_vv*abs(v);
    d_r = -N_r - N_rr*abs(r);

    %%%%%%%%% Dynamics
    ceq = [
        DiffMat*x - u.*cos(yaw) + v.*sin(yaw) - Vcx;
        DiffMat*y - u.*sin(yaw) - v.*cos(yaw) - Vcy;
        DiffMat*yaw - r;
        DiffMat*u - 1/m_u*(tau_u + m_v*v.*r - d_u.*u+fu);
        DiffMat*v - 1/m_v*(-m_u*u.*r - d_v.*v+fv);
        DiffMat*r - 1/m_r*(tau_r + m_uv*u.*v - d_r.*r+fr);
    ];
%     ceq = [
%         DiffMat*x - u.*cos(yaw) + v.*sin(yaw) - Vcx;
%         DiffMat*y - u.*sin(yaw) - v.*cos(yaw) - Vcy;
%         DiffMat*yaw - r;
%         DiffMat*u + 1/m_u*(-m_v.*v.*r+d_u-tau_u-fu);
%         DiffMat*v + 1/m_v*(m_u.*u.*r+d_v.*v-fv);
%         DiffMat*r + 1/m_r*(-m_uv.*u.*v+d_r.*r-tau_r-fr); 
%     ];
    
    c = [];
    
end

function J = costfun_single(X,constants) 
%     J = sum(X(:,7).^2)+sum(X(:,8).^2);
    %J = sum(X(:,7).^2)/constants.inputbounds(2,1) + ...
    %sum(X(:,7).^2)/constants.inputbounds(2,2);
    J=sum(BernsteinMul(X(:,4),X(:,7)).^2)+sum(BernsteinMul(X(:,6),X(:,8)).^2);
    % another cost to implement: integral of (force*vel + torque*angvel)
end

function xinit = init_guess(constants)

    N = constants.N; 

    xinit = zeros((constants.N-1)*constants.numvars,constants.Nv);
    for i = 1:constants.Nv
        x = linspace(constants.xi(i,1),constants.xf(i,1),N-1).';
        y = linspace(constants.xi(i,2),constants.xf(i,2),N-1).';
        yaw = ones(N-1,1).*atan2(y(end)-y(1),x(end)-x(1));
        u = zeros(N-1,1);
        v = zeros(N-1,1);
        r = zeros(N-1,1);
        xinit(:,i) = [x;y;yaw;u;v;r];
    end
    %xinit = xinit(:);
    xinit = [xinit ; zeros((constants.N+1)*constants.numinputs,constants.Nv)];

end
