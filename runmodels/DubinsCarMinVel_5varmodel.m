clear all; %#ok<CLALL>

addpath('..\Bernstein');
addpath('..\BeBOT_lib');
addpath('..\TrajecOptimLib');

% parameters

% 1 vehicle no constraints

% constants.T = 10;
% constants.xi = [ 0 0 0 1 0];
% constants.xf = [ 5 5 pi/2 1 0];
% % constants.obstacles_circles = [5,0,3];
% constants.plotboatsize = 0.7;
% constants.N = 15;
% constants.uselogbar = true; % use log barrier func


% constants.T = 20; % time
% constants.xi = [5 3 0 1 0]; % x y psi v w
% constants.xf = [0 0 0 1 0]; % x y psi v w
% constants.obstacles = [];
% constants.obstacles_circles = [];% = surroundhull(constants.obstacles);

% constants.T = 30;
% constants.xi  = [-10 40 0 1.1 0]; % x y psi v w
% constants.xf = [0 0 0 1.2 0]; % x y psi v w
% constants.obstacles = [];
% constants.obstacles_circles = [];
% constants.N = 100;
% constants.plotboatsize=2;

% 1 vehicle 1 obstacle
% constants.T = 15;
% constants.xi = [5 3 0 1 0];
% constants.xf = [0 0 0 1 0];
% constants.obstacles = [  1 1 ; 1 2 ; 2 2 ; 2 1 ];
% constants.obstacles_circles = [ 1.5, 1.5, 0.5];
% constants.plotboatsize = 0.5;

constants.T = 10;
constants.xi = [ -5 0 0 1 0]; % x y psi v r
constants.xf = [  5 0 0 1 0];
constants.N = 30;
% constants.obstacles = [ -0.5 -0.5; -0.5 0.5; 0.5 0.5; 0.5 -0.5 ];
constants.obstacles_circles = [0, 0, 1];
% constants.uselogbar = true;
% constants.useeqlogbar = true;

% 2 vehicles no obstacles
% constants.T = 15;
% constants.xi = [
%     0 5 0 1 0
%     5 0 pi/2 1 0
% ];
% constants.xf = [
%     10 5 0 1 0
%     5 10 pi/2 1 0 
% ];

% 3 vehicles 1 circle obstacle
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
% constants.obstacles_circles = [ 0 0 3]; % x y r

% common parameters
constants.min_dist_int_veh = .9;

constants.statebounds = [
    -Inf, -Inf, -2*pi, 0, -pi/4;
    Inf, Inf, 2*pi, 5, pi/4;
];

% functions
constants.costfun_single = @costfun_single;
constants.dynamics = @dynamicsDubin;
constants.init_guess = @init_guess;
constants.recoverxy = @recoverplot;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [xOut,JOut] = run_problem(constants);
[xOut,JOut] = run_problem_progressive_n(constants);
constants = processconstants(constants);
disp(['The final cost is ', num2str(JOut)])

%%

%constraint_evaluation(xOut,constants);
plot_xy(xOut,constants);
% plot some boats
points = BernsteinEval(xOut,constants.T,linspace(0,constants.T,10));
for v = 1:size(constants.xi, 1)
    for i = 1:10
        plotboat(points(i,1,v),points(i,2,v),points(i,3,v),constants.plotboatsize);
    end
end

%%
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

    c = [];

end


function J = costfun_single(X,constants)
    persistent nchoosek_mat N
    if isempty(nchoosek_mat) || constants.N ~= N
        N = constants.N;
        nchoosek_mat = nchoosek_mod_mat(2*N+3);
    end
    v = X(:,4);
    w = X(:,5);
    dv = constants.DiffMat*v;
    dw = constants.DiffMat*w;
%     J = sum(dv.^2)+2*sum(dw.^2);
    J = BernsteinIntegr(BernsteinMul(dv, dv, nchoosek_mat)+2*BernsteinMul(dw, dw, nchoosek_mat), constants.T);
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
