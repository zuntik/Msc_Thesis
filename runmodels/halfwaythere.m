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
CONSTANTS.T = 15; % time
CONSTANTS.xi = [0 0 0    1 0]; % x y yaw u r
%CONSTANTS.xf = [5 5 pi/2 1 0]; % x y yaw u r
% CONSTANTS.xf = [5 2 pi/6 1 0];
CONSTANTS.xf = [ 15 0 0 1 0 ];
CONSTANTS.N = 30; % order
CONSTANTS.obstacles = [];
CONSTANTS.obstacles_circles = [];% = surroundhull(CONSTANTS.obstacles);

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
CONSTANTS.dynamics = @dynamics;
CONSTANTS.init_guess = @rand_init_guess;
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

    tau_u = @(t) BernsteinEval(X(:,6),CONSTANTS.T,t);
    tau_r = @(t) BernsteinEval(X(:,7),CONSTANTS.T,t);
    
    [t,xy] = ode45(@(t,xy)odefunc(t,xy,tau_u,tau_r,CONSTANTS.MODELPARAMS), [0 CONSTANTS.T], X(1,1:5));

    function dydt = odefunc(t,y,tau_u,tau_r,MODELPARAMS)

        % load the parameters
        %%%%%% Coefficients
        X_u = MODELPARAMS.X_u;
        N_r = MODELPARAMS.N_r;
        X_uu = MODELPARAMS.X_uu;
        N_rr = MODELPARAMS.N_rr;

        %%%%%% masses
        m_u = MODELPARAMS.m_u;
        m_r = MODELPARAMS.m_r;

        yaw = y(3);
        u = y(4);
        r = y(5);

        %%%%%%%%% drag
        d_u = -X_u - X_uu*abs(u);
        d_r = -N_r - N_rr*abs(r);

        dydt = zeros(5,1);
        dydt(1) = u*cos(yaw); %dx
        dydt(2) = u*sin(yaw); % dy
        dydt(3) = r; % dyaw
        dydt(4) = 1/m_u*(tau_u(t) - d_u*u); %du
        dydt(5) = 1/m_r*(tau_r(t) - d_r*r); % dr  
    end

end

function [c,ceq] = dynamics(X,CONSTANTS)

    DiffMat = CONSTANTS.DiffMat;
    
    % states
    x = X(:,1);
    y = X(:,2);
    yaw = X(:,3);
    u = X(:,4);
    r = X(:,5);
    % inputs
    tau_u = X(:,6);
    tau_r = X(:,7); 

    % load the parameters
    %%%%%% Coefficients
%     mass = CONSTANTS.MODELPARAMS.mass;
%     I_z = CONSTANTS.MODELPARAMS.I_z;
%     X_dot_u = CONSTANTS.MODELPARAMS.X_dot_u;
%     Y_dot_v = CONSTANTS.MODELPARAMS.Y_dot_v;
%     N_dot_r = CONSTANTS.MODELPARAMS.N_dot_r;
    X_u = CONSTANTS.MODELPARAMS.X_u;
    N_r = CONSTANTS.MODELPARAMS.N_r;
    X_uu = CONSTANTS.MODELPARAMS.X_uu;
    N_rr = CONSTANTS.MODELPARAMS.N_rr;

    %%%%%% masses
    m_u = CONSTANTS.MODELPARAMS.m_u;
    m_r = CONSTANTS.MODELPARAMS.m_r;

    %%%%%%%%% drag
    d_u = -X_u - X_uu*abs(u);
    d_r = -N_r - N_rr*abs(r);

    %%%%%%%%% Dynamics
    ceq = [
        DiffMat*x - u.*cos(yaw);
        DiffMat*y - u.*sin(yaw);
        DiffMat*yaw - r;
        DiffMat*u - 1/m_u*(tau_u - d_u.*u);
        DiffMat*r - 1/m_r*(tau_r - d_r.*r);
    ];   

    c = [];
    
end

function J = costfun_single(X,CONSTANTS) 
    J = sum(X(:,6).^2)+sum(X(:,7).^2);
end

function xinit = rand_init_guess(CONSTANTS) %#ok<*DEFNU>
    xinit = rand((CONSTANTS.numvars*(CONSTANTS.N-1)+...
        CONSTANTS.numinputs*(CONSTANTS.N+1))*CONSTANTS.Nv,1);
end
