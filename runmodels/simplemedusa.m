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

CONSTANTS.MODELPARAMS=MODELPARAMS;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% runtime parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % 1 vehicle no constraints
% CONSTANTS.T = 15; % time
% CONSTANTS.xi = [0 0 0    1 0 0]; % x y yaw u v r
% %CONSTANTS.xf = [5 5 pi/2 1 0 0]; % x y yaw u v r
% CONSTANTS.xf = [5 2 pi/6 1 0 0];
% CONSTANTS.N = 50; % order
% CONSTANTS.obstacles = [];
% CONSTANTS.obstacles_circles = [];% = surroundhull(CONSTANTS.obstacles);

% 3 vehicles 1 circle obstacle
CONSTANTS.N = 40;
CONSTANTS.T = 20;
CONSTANTS.xi = [
    -10 4 0 1 0 0
%     -10 -4 0 1 0 0
%     -10 0 0 1 0 0
];
CONSTANTS.xf = [
    10 -1 0 1 0 0
%     10 1 0 1 0 0
%     10 0 0 1 0 0
];
CONSTANTS.T = 15;
CONSTANTS.obstacles = [];
CONSTANTS.obstacles_circles = [ 0 0 3]; % x y r


% common parameters
CONSTANTS.min_dist_intervehicles = 3;
CONSTANTS.DiffMat = BernsteinDerivElevMat(CONSTANTS.N,CONSTANTS.T);
CONSTANTS.EvalMat = BernsteinCtrlPntsEval(CONSTANTS.N);
CONSTANTS.BigElevMat = BernsteinDegrElevMat(CONSTANTS.N,CONSTANTS.N*10);
CONSTANTS.numvars = size(CONSTANTS.xi,2);
CONSTANTS.numinputs = 2;
CONSTANTS.Nv = size(CONSTANTS.xi,1);%number of vehicles
CONSTANTS.uselogbar = false; % use log barrier func
CONSTANTS.usesigma = true; % a variable for the usage of log barrier func

% functions
CONSTANTS.costfun_single = @costfun_single;
CONSTANTS.dynamics = @dynamicsmedusa;
CONSTANTS.init_guess = @init_guess;
CONSTANTS.recoverxy = @recoverplot;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xOut = run_problem(CONSTANTS);
%%
constraint_evaluation(xOut,CONSTANTS);
plot_xy(xOut,CONSTANTS);

times = linspace(0,CONSTANTS.T,10);

points = BernsteinEval(xOut,CONSTANTS.T,times);
for i = 1:10
    plotboat(points(i,1),points(i,2),points(i,3),0.5);
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
    BernsteinPlot(sum((xOut(:,1:2)-CONSTANTS.obstacles_circles(1:2)).^2,2),CONSTANTS.T);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t,xy] = recoverplot(X,CONSTANTS)

    tau_u = @(t) BernsteinEval(X(:,7),CONSTANTS.T,t);
    tau_r = @(t) BernsteinEval(X(:,8),CONSTANTS.T,t);
    
    [t,xy] = ode45(@(t,xy)odefunc(t,xy,tau_u,tau_r,CONSTANTS.MODELPARAMS), [0 CONSTANTS.T], X(1,1:6));

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

function [c,ceq] = dynamicsmedusa(X,CONSTANTS)

    DiffMat = CONSTANTS.DiffMat;
    
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
%     mass = CONSTANTS.MODELPARAMS.mass;
%     I_z = CONSTANTS.MODELPARAMS.I_z;
%     X_dot_u = CONSTANTS.MODELPARAMS.X_dot_u;
%     Y_dot_v = CONSTANTS.MODELPARAMS.Y_dot_v;
%     N_dot_r = CONSTANTS.MODELPARAMS.N_dot_r;
    X_u = CONSTANTS.MODELPARAMS.X_u;
    Y_v = CONSTANTS.MODELPARAMS.Y_v;
    N_r = CONSTANTS.MODELPARAMS.N_r;
    X_uu = CONSTANTS.MODELPARAMS.X_uu;
    Y_vv = CONSTANTS.MODELPARAMS.Y_vv;
    N_rr = CONSTANTS.MODELPARAMS.N_rr;

    %%%%%% masses
    m_u = CONSTANTS.MODELPARAMS.m_u;
    m_v = CONSTANTS.MODELPARAMS.m_v;
    m_r = CONSTANTS.MODELPARAMS.m_r;
    m_uv = CONSTANTS.MODELPARAMS.m_uv;

    %%%%%% Constants
    fu = CONSTANTS.MODELPARAMS.fu;
    fv = CONSTANTS.MODELPARAMS.fv;
    fr = CONSTANTS.MODELPARAMS.fr;
    Vcx = CONSTANTS.MODELPARAMS.Vcx;
    Vcy = CONSTANTS.MODELPARAMS.Vcy;

    %%%%%%%%% drag
    d_u = -X_u - X_uu*abs(u);
    d_v = -Y_v - Y_vv*abs(v);
    d_r = -N_r - N_rr*abs(r);

    %%%%%%%%% Dynamics
    ceq = [
        DiffMat*x - u.*cos(yaw) - v.*sin(yaw) + Vcx;
        DiffMat*y - u.*sin(yaw) + v.*cos(yaw) + Vcy;
        DiffMat*yaw - r;
        DiffMat*u - 1/m_u*(tau_u + m_v*v.*r - d_u.*u+fu);
        DiffMat*v - 1/m_v*(-m_u*u.*r - d_v.*v+fv);
        DiffMat*r - 1/m_r*(tau_r + m_uv*u.*v - d_r.*r+fr);
    ];

    c = [];
    
end

function J = costfun_single(X,CONSTANTS) 
%     v = X(:,4);
%     w = X(:,5);
%     a = CONSTANTS.DiffMat*v;
%     J = sum(a.^2)+2*sum(w.^2);
    J = sum(X(:,7).^2)+2*sum(X(:,6).^2);
end

function xinit = init_guess(CONSTANTS)

    N = CONSTANTS.N; 

    xinit = zeros((CONSTANTS.N-1)*CONSTANTS.numvars,CONSTANTS.Nv);
    for i = 1:CONSTANTS.Nv
        x = linspace(CONSTANTS.xi(i,1),CONSTANTS.xf(i,1),N-1).';
        y = linspace(CONSTANTS.xi(i,2),CONSTANTS.xf(i,2),N-1).';
        yaw = ones(N-1,1).*atan2(y(end)-y(1),x(end)-x(1));
        u = ones(N-1,1);
        v = zeros(N-1,1);
        r = zeros(N-1,1);
        xinit(:,i) = [x;y;yaw;u;v;r];
    end
    %xinit = xinit(:);
    xinit = [xinit ; rand((CONSTANTS.N+1)*CONSTANTS.numinputs,CONSTANTS.Nv)];

end

function xinit = rand_init_guess(CONSTANTS) %#ok<*DEFNU>
    xinit = rand((CONSTANTS.numvars*(CONSTANTS.N-1)+...
        CONSTANS.numinputs*(CONSTANTS.N+1))*CONSTANTS.Nv,1);
end


