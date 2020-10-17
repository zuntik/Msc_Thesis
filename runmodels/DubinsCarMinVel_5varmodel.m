clear all; %#ok<CLALL>

addpath('..\Bernstein');
addpath('..\BeBOT_lib');
addpath('..\TrajecOptimLib');

% parameters

% 1 vehicle no constraints

constants.T = 10;
constants.xi = [ 0 0 0 1 0];
constants.xf = [ 5 5 pi/2 1 0];
constants.N = 15;
constants.obstacles_circles = [5,0,3];

% constants.T = 15; % time
% constants.xi = [5 3 0 1 0]; % x y psi v w
% constants.xf = [0 0 0 1 0]; % x y psi v w
% constants.N = 50; % order
% constants.obstacles = [];
% constants.obstacles_circles = [];% = surroundhull(constants.obstacles);

% constants.T = 30;
% constants.xi  = [-10 40 0 1.1 0]; % x y psi v w
% constants.xf = [0 0 0 1.2 0]; % x y psi v w
% constants.N = 50; % order
% constants.obstacles = [];
% constants.obstacles_circles = [];

% 1 vehicle 1 obstacle
% constants.T = 15;
% constants.xi = [5 3 0 1 0];
% constants.xf = [0 0 0 1 0];
% constants.obstacles = [  1 1 ; 1 2 ; 2 2 ; 2 1 ] + [100 100];
% constants.obstacles_circles = [];
% constants.N = 13;

% constants.T = 10;
% constants.xi = [ -5 0 0 1 0];
% constants.xf = [  5 0 0 1 0];
% constants.N = 13;
% constants.obstacles = [ -0.5 -0.5; -0.5 0.5; 0.5 0.5; 0.5 -0.5 ];
% constants.obstacles_circles = [];

% 2 vehicles no obstacles
% constants.N = 30;
% constants.T = 15;
% constants.xi = [
%     0 5 0 1 0
%     5 0 pi/2 1 0
% ];
% constants.xf = [
%     10 5 0 1 0
%     5 10 pi/2 1 0 
% ];
% constants.obstacles = [];
% constants.obstacles_circles = [];

% 3 vehicles 1 circle obstacle
% constants.N = 40;
% constants.T = 15;
% constants.xi = [
%     -10 4 0 1 0
%     -10 -4 0 1 0
%     -10 0 0 1 0
% ];
% constants.xf = [
%     10 -1 0 1 0
%     10 1 0 1 0
%     10 0 0 1 0
% ];
% constants.obstacles = [];
% constants.obstacles_circles = [ 0 0 3]; % x y r

% common parameters
constants.min_dist_int_veh = .9;
% constants.numvars = size(constants.xi,2);
constants.Nv = size(constants.xi,1);%number of vehicles
% constants.uselogbar = true; % use log barrier func
% constants.usesigma = true; % a variable for the usage of log barrier func

% functions
constants.costfun_single = @costfun_single;
constants.dynamics = @dynamicsDubin;
% constants.init_guess = @init_guess;
constants.recoverxy = @recoverplot;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[xOut,JOut] = run_problem(constants);

disp(['The final cost is ', num2str(JOut)])
%%
%constraint_evaluation(xOut,constants);
plot_xy(xOut,constants);

if constants.Nv == 2
    figure
    X_diff = xOut(:,1:2,1)-xOut(:,1:2,2);
    X_diff_2 = BernsteinPow(X_diff,2);
    Mag = sqrt(sum(X_diff_2,2));
    BernsteinPlot(Mag,constants.T);
end

if (isfield(constants,'obstacles_circles') && ~isempty(constants.obstacles_circles))...
        && size(constants.obstacles_circles,3) == 1 && constants.Nv == 1
    figure, grid on
    BernsteinPlot(sum((xOut(:,1:2)-constants.obstacles_circles(1:2)).^2,2),constants.T);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t,xy] = recoverplot(X,constants)

    v = @(t) BernsteinEval(X(:,4),constants.T,t);
    w = @(t) BernsteinEval(X(:,5),constants.T,t);
    
    [t,xy] = ode45(@(t,xy)odefunc(t,xy,v,w), [0 constants.T], X(1,1:3));

    function dydt = odefunc(t,y,v,w)

        dydt = zeros(3,1);
        dydt(1) = v(t)*cos(y(3));%x
        dydt(2) = v(t)*sin(y(3));%y
        dydt(3) = w(t);%psi
        
    end

end


function [c,ceq] = dynamicsDubin(X,constants)

    DiffMat = constants.DiffMat;
    
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
    %c = constants.maxtorque- Xout(:,8);
    
end


function J = costfun_single(X,constants)
    v = X(:,4);
    r = X(:,5);
    a = constants.DiffMat*v;
    J = sum(a.^2)+2*sum(r.^2);
end


function xinit = init_guess(constants)

    N = constants.N; 

    xinit = zeros((constants.N-1)*constants.numvars,constants.Nv);
    for i = 1:constants.Nv
        x = linspace(constants.xi(i,1),constants.xf(i,1),N-1).';
        y = linspace(constants.xi(i,2),constants.xf(i,2),N-1).';
        %psi = -ones(N-1,1);
        psi = ones(N-1,1).*atan2(y(end)-y(1),x(end)-x(1));
        v = ones(N-1,1);
        w = zeros(N-1,1);
        xinit(:,i) = [x;y;psi;v;w];
    end
    %xinit = xinit(:);
    
end


function xinit = rand_init_guess(constants) %#ok<*DEFNU>

    N = constants.N; 

    x = rand(N-1,1);
    y = rand(N-1,1);
    psi = rand(N-1,1);
    v = rand(N-1,1);
    w = rand(N-1,1);

    xinit = [x;y;psi;v;w];
    xinit = xinit(:);
    
end
