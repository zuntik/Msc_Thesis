

%dubin

%%%%%%%%%%%%%%%% Example 1 %%%%%%%%%%%%%%%%
% 1 vehicle

% constants.T = 10;
% constants.xi = [ 0 0 0 1 0];
% constants.xf = [ 5 5 pi/2 1 0];
% % constants.obstacles_circles = [5,0,3];
% constants.obstacles = [ 3 2; 3 0; 6 0; 6 2];
% constants.N = 15;
% % constants.uselogbar = true;
% constants.useeqlogbar = false;

%%%%%%%%%%%%%%%% Example 2 %%%%%%%%%%%%%%%%
% 1 vehicle

% constants.T = 20; % time
% constants.xi = [5 3 0 1 0]; % x y psi v w
% constants.xf = [0 0 0 1 0]; % x y psi v w
% constants.obstacles = [  1 1 ; 1 2 ; 2 2 ; 2 1 ];
% % constants.obstacles_circles = [ 1.5, 1.5, 0.5];
% constants.N = 20;

%%%%%%%%%%%%%%%% Example 3 %%%%%%%%%%%%%%%%
% 1 vehicle

% constants.T = 10;
% constants.xi = [ -5 0 0 1 0]; % x y psi v r
% constants.xf = [  5 0 0 1 0];
% constants.N = 10;
% % constants.obstacles = 0.5 *[ -1 -1; -1 1; 1 1; 1 -1 ];
% constants.obstacles_circles = [0, 0, 1];
% constants.uselogbar = true;
% constants.useeqlogbar = false;



%%%%%%%%%%%%%%%% Example 4 %%%%%%%%%%%%%%%%
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


%%%%%%%%%%%%%%%% Example 6  - Many vehicles!! %%%%%%%%%%%%%%%%
% constants.T = 100;
% constants.N = 10;
% constants.Nv = 10;
% constants.xi = repmat([0  0 0 0.1 0] ,[constants.Nv 1]);
% constants.xf = repmat([14 0 0 0.1 0] ,[constants.Nv 1]);
% constants.xi(:, 2) = 2*(0:1:constants.Nv -1).';
% constants.xf(:, 2) = constants.xi(randperm(size(constants.xi,1)),2);

%%%%%%%%%%%%%%%% State variable bounds %%%%%%%%%%%%%%%%
% constants.statebounds = [
%     -Inf, -Inf, -2*pi, 0, -pi/4; % inferior bounds of each state variable
%     Inf, Inf, 2*pi, 5, pi/4; % superior bounds of each state variable
% ];


% medusa

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
