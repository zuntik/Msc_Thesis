clear all; %#ok<CLALL>

addpath('..\Bernstein');
addpath('..\BeBOT_lib');
addpath('..\TrajecOptimLib');

% parameters

% 1 vehicle no constraints

CONSTANTS.T = 10;
CONSTANTS.xi = [ 0 0 0 1 0];
CONSTANTS.xf = [ 5 5 pi/2 1 0];
CONSTANTS.N = 15;
CONSTANTS.obstacles = [];
CONSTANTS.obstacles_circles = [];
%CONSTANTS.obstacles_circles = [5,0,3];

% CONSTANTS.T = 15; % time
% CONSTANTS.xi = [5 3 0 1 0]; % x y psi v w
% CONSTANTS.xf = [0 0 0 1 0]; % x y psi v w
% CONSTANTS.N = 50; % order
% CONSTANTS.obstacles = [];
% CONSTANTS.obstacles_circles = [];% = surroundhull(CONSTANTS.obstacles);

% CONSTANTS.T = 30;
% CONSTANTS.xi  = [-10 40 0 1.1 0]; % x y psi v w
% CONSTANTS.xf = [0 0 0 1.2 0]; % x y psi v w
% CONSTANTS.N = 50; % order
% CONSTANTS.obstacles = [];
% CONSTANTS.obstacles_circles = [];

% 1 vehicle 1 obstacle
% CONSTANTS.T = 15;
% CONSTANTS.xi = [5 3 0 1 0];
% CONSTANTS.xf = [0 0 0 1 0];
% CONSTANTS.obstacles = [  1 1 ; 1 2 ; 2 2 ; 2 1 ] + [100 100];
% CONSTANTS.obstacles_circles = [];
% CONSTANTS.N = 13;

% CONSTANTS.T = 10;
% CONSTANTS.xi = [ -5 0 0 1 0];
% CONSTANTS.xf = [  5 0 0 1 0];
% CONSTANTS.N = 13;
% CONSTANTS.obstacles = [ -0.5 -0.5; -0.5 0.5; 0.5 0.5; 0.5 -0.5 ];
% CONSTANTS.obstacles_circles = [];

% 2 vehicles no obstacles
% CONSTANTS.N = 30;
% CONSTANTS.T = 15;
% CONSTANTS.xi = [
%     0 5 0 1 0
%     5 0 pi/2 1 0
% ];
% CONSTANTS.xf = [
%     10 5 0 1 0
%     5 10 pi/2 1 0 
% ];
% CONSTANTS.obstacles = [];
% CONSTANTS.obstacles_circles = [];

% 3 vehicles 1 circle obstacle
% CONSTANTS.N = 40;
% CONSTANTS.T = 15;
% CONSTANTS.xi = [
%     -10 4 0 1 0
%     -10 -4 0 1 0
%     -10 0 0 1 0
% ];
% CONSTANTS.xf = [
%     10 -1 0 1 0
%     10 1 0 1 0
%     10 0 0 1 0
% ];
% CONSTANTS.obstacles = [];
% CONSTANTS.obstacles_circles = [ 0 0 3]; % x y r

% common parameters
CONSTANTS.min_dist_intervehicles = .9;
CONSTANTS.DiffMat = BernsteinDerivElevMat(CONSTANTS.N,CONSTANTS.T);
CONSTANTS.EvalMat = BernsteinCtrlPntsEval(CONSTANTS.N);
CONSTANTS.BigElevMat = BernsteinDegrElevMat(CONSTANTS.N,CONSTANTS.N*10);
CONSTANTS.numvars = size(CONSTANTS.xi,2);
CONSTANTS.numinputs = 0;
CONSTANTS.Nv = size(CONSTANTS.xi,1);%number of vehicles
CONSTANTS.uselogbar = true; % use log barrier func
CONSTANTS.usesigma = true; % a variable for the usage of log barrier func

% functions
CONSTANTS.costfun_single = @costfun_single;
CONSTANTS.dynamics = @dynamics5vars;
CONSTANTS.init_guess = @init_guess;
CONSTANTS.recoverxy = @recoverplot;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[xOut,JOut] = run_problem(CONSTANTS);

disp(['The final cost is ', num2str(JOut)])
%%
%constraint_evaluation(xOut,CONSTANTS);
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
function [t,xy] = recoverplot(X,CONSTANTS)

    v = @(t) BernsteinEval(X(:,4),CONSTANTS.T,t);
    w = @(t) BernsteinEval(X(:,5),CONSTANTS.T,t);
    
    [t,xy] = ode45(@(t,xy)odefunc(t,xy,v,w), [0 CONSTANTS.T], X(1,1:3));

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
    %c = [-v+0.2;psi-pi;-psi-pi]; % this one is venanzio's
    %c = [ -v; psi-pi-pi/3; -psi-pi-pi/3];
    %c = -v + 0.2 ;
    c = [];
    %c = CONSTANTS.maxtorque- Xout(:,8);
    
end


function J = costfun_single(X,CONSTANTS)
    v = X(:,4);
    w = X(:,5);
    a = CONSTANTS.DiffMat*v;
    J = sum(a.^2)+2*sum(w.^2);
end


function xinit = init_guess(CONSTANTS)

    N = CONSTANTS.N; 

    xinit = zeros((CONSTANTS.N-1)*CONSTANTS.numvars,CONSTANTS.Nv);
    for i = 1:CONSTANTS.Nv
        x = linspace(CONSTANTS.xi(i,1),CONSTANTS.xf(i,1),N-1).';
        y = linspace(CONSTANTS.xi(i,2),CONSTANTS.xf(i,2),N-1).';
        %psi = -ones(N-1,1);
        psi = ones(N-1,1).*atan2(y(end)-y(1),x(end)-x(1));
        v = ones(N-1,1);
        w = zeros(N-1,1);
        xinit(:,i) = [x;y;psi;v;w];
    end
    %xinit = xinit(:);
    
end


function xinit = rand_init_guess(CONSTANTS) %#ok<*DEFNU>

    N = CONSTANTS.N; 

    x = rand(N-1,1);
    y = rand(N-1,1);
    psi = rand(N-1,1);
    v = rand(N-1,1);
    w = rand(N-1,1);

    xinit = [x;y;psi;v;w];
    xinit = xinit(:);
    
end
