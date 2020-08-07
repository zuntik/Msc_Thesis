clear all; %#ok<CLALL>

addpath('..\Bernstein');
addpath('..\BeBOT_lib');
addpath('..\TrajecOptimLib');

% parameters

% 1 vehicle no constraints
CONSTANTS.T = 15; % time
CONSTANTS.xi = [5 3 0 1 0]; % x y psi v w
CONSTANTS.xf = [0 0 0 1 0]; % x y psi v w
CONSTANTS.N = 50; % order
CONSTANTS.obstacles = [];
CONSTANTS.obstacles_circles = [];% = surroundhull(CONSTANTS.obstacles);


% common parameters
CONSTANTS.min_dist_intervehicles = 3;
CONSTANTS.DiffMat = BernsteinDerivElevMat(CONSTANTS.N,CONSTANTS.T);
CONSTANTS.EvalMat = BernsteinCtrlPntsEval(CONSTANTS.N);
CONSTANTS.BigElevMat = BernsteinDegrElevMat(CONSTANTS.N,CONSTANTS.N*10);
CONSTANTS.numvars = size(CONSTANTS.xi,2);
CONSTANTS.numinpus = 2;
CONSTANTS.Nv = size(CONSTANTS.xi,1);%number of vehicles
CONSTANTS.uselogbar = false; % use log barrier func
CONSTANTS.usesigma = true; % a variable for the usage of log barrier func

% functions
CONSTANTS.costfun_single = @costfun_single;
CONSTANTS.dynamics = @dynamics5vars;
CONSTANTS.init_guess = @init_guess;
CONSTANTS.recoverxy = @recoverplot;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xOut = run_problem(CONSTANTS);
constraint_evaluation(xOut,CONSTANTS);
plot_xy(xOut,CONSTANTS);

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

    %%%%%% Constants
    fu  = 0;
    fv  = 0;
    fr  = 0;
    Vcx = 0;
    Vcy = 0;
    
    % %%%%%% Coefficients
    % mass = 17.0;
    % I_z = 4.14;
    % X_dot_u = -20;
    % Y_dot_v = -1.3175;
    % N_dot_r =  -8.69;
    % X_u = -0.2;
    % Y_v = -55.117;
    % N_r = -4.14;
    % X_uu = -25; 
    % Y_vv = -101.2776;
    % N_rr = -6.23;

    %%%%%% Coefficients
    mass = 17.0;
    I_z = 1;
    X_dot_u = -20;
    Y_dot_v = -30;%-1.3175;
    N_dot_r = -8.69;% -0.5;
    X_u = -0.2;
    Y_v = -50;
    N_r = -4.14; %-0.1
    X_uu = -25; 
    Y_vv = -0.01;%-101.2776;
    N_rr = -6.23; %-21

    %%%%%% masses
    m_u = mass - X_dot_u;
    m_v = mass - Y_dot_v;
    m_r = I_z - N_dot_r;
    m_uv = m_u - m_v;

    %%%%%%%%% drag

    d_u = -X_u - X_uu*abs(u);
    d_v = -Y_v - Y_vv*abs(v);
    d_r = -N_r - N_rr*abs(r);

    %%%%%%%%% Dynamics
    ceq = [
        DiffMat*u - 1/m_u*(tau_u + m_v*v*r - d_u*u+fu);
        DiffMat*v - 1/m_v*(-m_u*u*r - d_v*v+fv);
        DiffMat*r - 1/m_r*(tau_r + m_uv*u*v - d_r*r+fr);
        DiffMat*x - u*cos(yaw) - v*sin(yaw) + Vcx;
        DiffMat*y - u*sin(yaw) + v*cos(yaw) + Vcy;
        DiffMat*yaw - r;
    ];

    c = [];
    
end

function J = costfun_single(X,CONSTANTS)
    v = X(:,4);
    w = X(:,5);
    a = CONSTANTS.DiffMat*v;
    J = sum(a.^2)+2*sum(w.^2);
    
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

