clear all; %#ok<CLALL>

addpath('..\Bernstein');
addpath('..\BeBOT_lib');
addpath('..\TrajecOptimLib');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% medusa model parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MODELPARAMS = struct();

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
constants.T = 15; % time
constants.xi = [0 0 0]; % x y yaw u v r
%constants.xf = [5 5 pi/2 1 0 0]; % x y yaw u v r
% constants.xf = [5 2 pi/6 1 0 0];
constants.xf = [ 0 6 pi ];
constants.N = 50; % order
constants.obstacles = [];
constants.obstacles_circles = [];% = surroundhull(constants.obstacles);
constants.ui = 1;
constants.uf = 1;
constants.vi = 0;
constants.vf = 0;
constants.ri = 0;
constants.rf = 0;


% common parameters
constants.min_dist_intervehicles = 3;
constants.numinputs = 0; % 2;

constants.statebounds = [
    % x   y    yaw
    -100 -100 -100 % lower state bounds
    100 100 100 % upper state bounds
    ];
constants.inputbounds = [
    % t_u t_r
    -100 -5; % lower input bounds
    100 5 % upper input bounds
    ];

% functions
constants.costfun_single = @costfun_single;
constants.dynamics = @dynamicsmedusa;
constants.init_guess = @init_guess;
constants.recoverxy = @recoverplot;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xOut = run_problem(constants);
%%
%constraint_evaluation(xOut,constants);
plot_xy(xOut,constants);

times = linspace(0,constants.T,10);

points = BernsteinEval(xOut,constants.T,times);
for i = 1:10
    plotboat(points(i,1),points(i,2),points(i,3),0.5);
end

if ~isempty(constants.obstacles_circles) && size(constants.obstacles_circles,3) == 1 && constants.Nv == 1
    figure, grid on
    BernsteinPlot(sum((xOut(:,1:2)-constants.obstacles_circles(1:2)).^2,2),constants.T);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t,xy] = recoverplot(X,constants)


%     tau_u = @(t)calc_tau_u(X, constants, t);
%     tau_r = @(t)calc_tau_r(X, constants, t);
    [~, ~, ~, tau_u, tau_r] = calcothers(X,constants);
    tau_u = @(t) BernsteinEval(tau_u,constants.T,t);
    tau_r = @(t) BernsteinEval(tau_r,constants.T,t);
    
    xi = [X(1,1:3), constants.ui, constants.vi, constants.ri];
    
    [t,xy] = ode45(@(t,xy)odefunc(t,xy,tau_u,tau_r,constants.MODELPARAMS), [0 constants.T], xi);

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

function [c,ceq] = dynamicsmedusa(X,constants)

    DiffMat = constants.DiffMat;
    
    % states
    x = X(:,1);
    y = X(:,2);
    yaw = X(:,3);
%     % inputs
%     tau_u = X(:,4);
%     tau_r = X(:,5); 

    % load the parameters
    %%%%%% Coefficients
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

    dx = DiffMat*x;
    dy = DiffMat*y;
    u = (dx-Vcx).*cos(yaw) + (dy-Vcy).*sin(yaw);
    v = -(dx-Vcx).*sin(yaw) + (dy-Vcy).*cos(yaw);
    r = DiffMat*yaw;

    %%%%%%%%% drag
    d_u = -X_u - X_uu*abs(u);
    d_v = -Y_v - Y_vv*abs(v);
    d_r = -N_r - N_rr*abs(r);
    
    tau_u = m_u*DiffMat*u - m_v*v.*r + d_u.*u;
    tau_r = m_r*DiffMat*r - m_uv * u.*v + d_r .*r;
    
    %%%%%%%%% Dynamics
    ceq = [
%         DiffMat*u - 1/m_u*(tau_u + m_v*v.*r - d_u.*u+fu);
        DiffMat*v - 1/m_v*(-m_u*u.*r - d_v.*v+fv);
%         DiffMat*r - 1/m_r*(tau_r + m_uv*u.*v - d_r.*r+fr);
    ];

    %%%%%%%%% Initial and final conditions
    ceq = [
        ceq;
        u(1)-constants.ui;
        u(end)-constants.uf;
        v(1)-constants.vi;
        v(end)-constants.vf;
        r(1)-constants.ri;
        r(end)-constants.rf;
    ];

    % limits to tau_u and tau_r
    c = [
        constants.inputbounds(1,1) - tau_u;
        tau_u - constants.inputbounds(2,1);
        constants.inputbounds(1,2) - tau_r;
        tau_r - constants.inputbounds(2,2);
    ];
    
end

function J = costfun_single(X,constants) 

    [~,~,~, tau_u, tau_r] = calcothers(X,constants);
    J = sum(tau_u.^2) + sum(tau_r.^2);

    %J = sum(X(:,4).^2)+sum(X(:,5).^2);
end

function xinit = init_guess(constants)

    N = constants.N; 

    xinit = zeros((constants.N-1)*constants.numvars,constants.Nv);
    for i = 1:constants.Nv
        x = linspace(constants.xi(i,1),constants.xf(i,1),N-1).';
        y = linspace(constants.xi(i,2),constants.xf(i,2),N-1).';
        yaw = ones(N-1,1).*atan2(y(end)-y(1),x(end)-x(1));
        xinit(:,i) = [x;y;yaw];
    end
    
    xinit = [xinit ; rand((constants.N+1)*constants.numinputs,constants.Nv)];

end


function [u,v,r,tau_u,tau_r] = calcothers(X,constants)

    DiffMat = constants.DiffMat;
    
    % states
    x = X(:,1);
    y = X(:,2);
    yaw = X(:,3);
%     % inputs
%     tau_u = X(:,4);
%     tau_r = X(:,5); 

    % load the parameters
    %%%%%% Coefficients
    X_u = constants.MODELPARAMS.X_u;
    N_r = constants.MODELPARAMS.N_r;
    X_uu = constants.MODELPARAMS.X_uu;
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

    dx = DiffMat*x;
    dy = DiffMat*y;
    u = (dx-Vcx).*cos(yaw) + (dy-Vcy).*sin(yaw);
    v = -(dx-Vcx).*sin(yaw) + (dy-Vcy).*cos(yaw);
    r = DiffMat*yaw;

    %%%%%%%%% drag
    d_u = -X_u - X_uu*abs(u);
    d_r = -N_r - N_rr*abs(r);
    
    tau_u = m_u*DiffMat*u - m_v*v.*r + d_u.*u;
    tau_r = m_r*DiffMat*r - m_uv * u.*v + d_r .*r;
end

function tau_u = calc_tau_u(X,constants, t)

    DiffMat = constants.DiffMat;
    
    % states
    x = X(:,1);
    y = X(:,2);
    yaw = X(:,3);

    % load the parameters
    %%%%%% Coefficients
    X_u = constants.MODELPARAMS.X_u;
    X_uu = constants.MODELPARAMS.X_uu;

    %%%%%% masses
    m_u = constants.MODELPARAMS.m_u;
    m_v = constants.MODELPARAMS.m_v;

    %%%%%% Constants
    Vcx = constants.MODELPARAMS.Vcx;
    Vcy = constants.MODELPARAMS.Vcy;

    dx = DiffMat*x;
    dy = DiffMat*y;
    dyaw = DiffMat*yaw;
    ddx = DiffMat*dx;
    ddy = DiffMat*dy;
    
    yaw = BernsteinEval(yaw, constants.T, t);
    dx  = BernsteinEval(dx,  constants.T, t);
    dy  = BernsteinEval(dy,  constants.T, t);
    dyaw = BernsteinEval(dyaw, constants.T, t);
    ddx = BernsteinEval(ddx, constants.T, t);
    ddy = BernsteinEval(ddy, constants.T, t);
    
    u =  (dx-Vcx).*cos(yaw) + (dy-Vcy).*sin(yaw);
    du = ddx.*cos(yaw)-(dx-Vcx).*sin(yaw).*dyaw+ddy.*sin(yaw)+(dy-Vcy).*cos(yaw).*dyaw;
    v = -(dx-Vcx).*sin(yaw) + (dy-Vcy).*cos(yaw);
    
    %%%%%%%%% drag
    d_u = -X_u - X_uu*abs(u);    
    tau_u = m_u*du - m_v*v.*dyaw + d_u.*u;

end

function tau_r = calc_tau_r(X, constants, t)

    DiffMat = constants.DiffMat;
    
    % states
    x = X(:,1);
    y = X(:,2);
    yaw = X(:,3);
    
    % load the parameters
    %%%%%% Coefficients
    N_r = constants.MODELPARAMS.N_r;
    N_rr = constants.MODELPARAMS.N_rr;

    %%%%%% masses
    m_r = constants.MODELPARAMS.m_r;
    m_uv = constants.MODELPARAMS.m_uv;

    %%%%%% Constants
    Vcx = constants.MODELPARAMS.Vcx;
    Vcy = constants.MODELPARAMS.Vcy;

    dx = DiffMat*x;
    dy = DiffMat*y;
    r = DiffMat*yaw;
    dr = DiffMat*r;
    
    yaw = BernsteinEval(yaw, constants.T, t);
    dx  = BernsteinEval(dx,  constants.T, t);
    dy  = BernsteinEval(dy,  constants.T, t);
    r   = BernsteinEval(r,   constants.T, t);
    dr  = BernsteinEval(dr, constants.T, t);
    
    u = (dx-Vcx).*cos(yaw) + (dy-Vcy).*sin(yaw);
    v = -(dx-Vcx).*sin(yaw) + (dy-Vcy).*cos(yaw);


    %%%%%%%%% drag
    d_r = -N_r - N_rr*abs(r);
    
    tau_r = m_r*dr - m_uv * u.*v + d_r .*r;
end
