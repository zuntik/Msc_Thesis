
% 1 vehicle no constraints
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
CONSTANTS.T = 15;
CONSTANTS.xi = [5 3 0 1 0];
CONSTANTS.xf = [0 0 0 1 0];
CONSTANTS.obstacles = [  1 1 ; 1 2 ; 2 2 ; 2 1 ] + [100 100];
CONSTANTS.obstacles_circles = [];
CONSTANTS.N = 13;


CONSTANTS.T = 10;
CONSTANTS.xi = [ -5 0 0 1 0];
CONSTANTS.xf = [  5 0 0 1 0];
CONSTANTS.N = 13;
CONSTANTS.obstacles = [ -0.5 -0.5; -0.5 0.5; 0.5 0.5; 0.5 -0.5 ];
CONSTANTS.obstacles_circles = [];

% 2 vehicles no obstacles
% CONSTANTS.T = 15;
% CONSTANTS.xi = [
%     0 5 0 1 0
%     5 0 pi/2 1 0
% ];
% CONSTANTS.xf = [
%     10 5 0 1 0
%     5 10 pi/2 1 0 
% ];

% 3 vehicles 1 circle obstacle
% CONSTANTS.T = 20;
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
% CONSTANTS.T = 15;
% CONSTANTS.obstacles = [];
% CONSTANTS.obstacles_circles = [ 0 0 3]; % x y r



CONSTANTS.min_dist_intervehicles = 3;
CONSTANTS.DiffMat = BernsteinDerivElevMat(CONSTANTS.N,CONSTANTS.T);
CONSTANTS.EvalMat = BernsteinCtrlPntsEval(CONSTANTS.N);
CONSTANTS.BigElevMat = BernsteinDegrElevMat(CONSTANTS.N,CONSTANTS.N*10);
CONSTANTS.numvars = size(CONSTANTS.xi,2);
CONSTANTS.Nv = size(CONSTANTS.xi,1);%number of vehicles
CONSTANTS.uselogbar = false; % use log barrier func
