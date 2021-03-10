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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% runtime parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fsolve(@func_for_circle, zeros(1, 4))
% disp('u, v, tau_u, tau_r')

% % 1 vehicle no constraints

constants.N = 20;

constants.T = 70; % time
constants.r = 10;
extraang = acos(0.8956/sqrt(0.0595^2+0.8956^2));
constants.xi = [constants.r, 0, pi/2+extraang, 0.8956, -0.0595, 2*pi/constants.T]; % x y yaw u v r
constants.xf = [constants.r, 0, pi/2+2*pi+extraang, 0.8956, -0.0595, 2*pi/constants.T ];
constants.desiredpoints=constants.r*[cos(linspace(0,2*pi,constants.N+1)).',...
    sin(linspace(0,2*pi,constants.N+1)).'];

% constants.T = 200;
% constants.desiredpoints=[ linspace(0,constants.T,constants.N+1).', 10.*sin(0.05*linspace(0,constants.T,constants.N+1)).'];
% constants.xi = [ 0, 0, 0.4630, 1.1177, 0, 0 ];
% constants.xf = [ 200, -5.4402, -0.4082, 1.0895, 0, 0 ];


% common parameters
constants.min_dist_intervehicles = 3;
constants.numinputs = 2;

constants.statebounds = [
    % x   y    yaw   u   v     r
    -1000 -1000 -1000 0 -1000 -.74 % lower state bounds
    1000 1000 1000 3 1000 .74 % upper state bounds
    ];
constants.inputbounds = [
    % t_u t_r
    0 -.113; % lower input bounds
    25.9 .113 % upper input bounds
    ];

% functions
constants.costfun_single = @costfun_single;
constants.dynamics = @dynamicsmedusa;
constants.init_guess = @exact_guess;
constants.recoverxy = @recoverplot; %_simulink;

%%%%%% Extras
constants.MODELPARAMS=MODELPARAMS;
constants.EvalMat = BernsteinEvalMat(constants.N,...
    constants.T, linspace(0, constants.T, 1000));

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xOut = run_problem(constants);
%%
constants = processconstants(constants);
plot_xy(xOut,constants);

figure
dxOut = BernsteinDeriv(xOut, constants.T);
speedangle = @(t) atan2(BernsteinEval(dxOut(:,2),constants.T,t), BernsteinEval(dxOut(:,1),constants.T,t));
sideslip = @(t) speedangle(t) - BernsteinEval(xOut(:,3),constants.T,t);
fplot(@(t)sideslip(t)/pi*180, [0, constants.T]);

constants.DiffMat = BernsteinDerivElevMat(constants.N,constants.T);
figure, hold on
tau_r_out = fplot(@(t)calc_tau_u(xOut,constants,t), [0,constants.T]);
BernsteinPlot(xOut(:,7),constants.T);
figure, hold on
tau_u_out = fplot(@(t)calc_tau_r(xOut,constants,t), [0,constants.T]);
BernsteinPlot(xOut(:,8),constants.T);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t,xy,tenpoints] = recoverplot(X,constants)

%     tau_u = @(t)calc_tau_u(X, constants, t);
%     tau_r = @(t)calc_tau_r(X, constants, t);
    tau_u = @(t) BernsteinEval(X(:,7),constants.T,t);
    tau_r = @(t) BernsteinEval(X(:,8),constants.T,t);

    [t,xy] = ode45(@(t,xy)odefunc(t,xy,tau_u,tau_r,constants.MODELPARAMS), [0 constants.T], X(1,1:6));

    tenpoints = BernsteinEval(X,constants.T,linspace(0,constants.T,10));

    function dydt = odefunc(t,y,tau_u,tau_r,MODELPARAMS)

        % load the parameters
        %%%%%% Coefficients
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
    u = X(:,4);
    v = X(:,5);
    r = X(:,6);
    % inputs
    tau_u = X(:,7);
    tau_r = X(:,8);

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
    c = [];

end

% function J = costfun_single(X,constants) 
%     t = linspace(0,2*pi,size(constants.EvalMat,1));
%     desiredpoints=constants.r*[cos(t).',sin(t).'];
%     J = sum((constants.EvalMat * X(:, 1:2) - desiredpoints).^2, 'all');
% end
function J = costfun_single(X,constants) 
    J = 0;
%     J = J + sum((X(:,3) - (pi/2+linspace(0,2*pi, constants.N+1)).').^2, 'all');
    J = J + sum((X(:, 1:2) - constants.desiredpoints(:,1:2)).^2, 'all');
end
    
function xinit = exact_guess(constants)
    xp = constants.desiredpoints(:,1);
    yp = constants.desiredpoints(:,2);
%     yaw = constants.desiredpoints(:,3);
%     xp = constants.r*cos(linspace(0,2*pi, constants.N+1)).';
%     yp = constants.r*sin(linspace(0,2*pi, constants.N+1)).';
    yaw = pi/2 + linspace(0, 2*pi, constants.N+1).';
    [u, v, r, tau_u, ~] = calcothers([xp, yp, yaw], constants);
    xvars = [xp, yp, yaw, u, v, r];
    xvars = xvars(2:end-1, :);
    tau_r = zeros(constants.N+1,1);
    xinput = [tau_u, tau_r];
    xinit = [xvars(:);xinput(:)];
end

function [u,v,r,tau_u,tau_r] = calcothers(X,constants)

    DiffMat = constants.DiffMat;

    % states
    x = X(:,1);
    y = X(:,2);
    yaw = X(:,3);

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

function d = xy_d(t)
    if t <= 10
        d = [ 0.9*t, 0 , 0];
    elseif t <= 5*pi/.9 + 10
%         d = [ 9, 10] + 10*[cos((t-10)*0.9/(5*pi)-pi/2), sin((t-10)*0.9/(5*pi)-pi/2)];
        d = [ 9, 10, 0] + 10*[cos((t-10)*0.09-pi/2), sin((t-10)*0.09-pi/2), (t-10)*0.09];
    elseif t <= 5*pi/.9 + 10 + 10
        d = [ 19, 10 + (t- (5*pi/.9 + 10)) * 0.9, pi/2];
    else
        error('should not be calculating outside this range');
    end
end



function F = func_for_circle(x)
% objective of fsolve is to make F = 0

T = 70;
radius = 10;
r = 2*pi/T;

mass = 17.0;
% I_z = 1;
X_dot_u = -20;
Y_dot_v = -30;%-1.3175;
% N_dot_r = -8.69;% -0.5;
X_u = -0.2;
Y_v = -50;
N_r = -4.14; %-0.1
X_uu = -25; 
Y_vv = -0.01;%-101.2776;
N_rr = -6.23; %-21

%%%%%% masses
m_u = mass - X_dot_u;
m_v = mass - Y_dot_v;
% m_r = I_z - N_dot_r;
m_uv = m_u - m_v;

u = x(1);
v = x(2);
tau_u = x(3);
tau_r = x(4);

d_u = -X_u -X_uu*abs(u);
d_v = -Y_v -Y_vv*abs(v);
d_r = -N_r -N_rr*abs(r);

F(1) = tau_u + m_v * v * r - d_u * u;
F(2) = m_u*u*r + d_v * v;
F(3) = tau_r + m_uv * u*v - d_r * r;
F(4) = u^2 + v^2 - radius^2*r^2;

end
