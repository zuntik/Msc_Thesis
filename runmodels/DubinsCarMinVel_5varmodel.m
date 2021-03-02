% clear all; %#ok<CLALL>

addpath('..\Bernstein');
addpath('..\BeBOT_lib');
addpath('..\TrajecOptimLib');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dubin car has five state variables: x,y,psi,u,v by that order
% they appear in xi, xf, and statebounds
% "constants" is what will define an optimisation probelm. it will have
% several predefined fields, however:
% to use log barrier function create uselogbar field and set to true
% to change order create N field and set to desired variable
constants = struct();

%%%%%%%%%%%%%%%% Example 0 %%%%%%%%%%%%%%%%
% constants.T = 80;
% %                x  y  psi  v  w
% constants.xi = [ 0  0  pi/4 .1 0 ]; % initial conds
% constants.xf = [ 10 10 pi/4 .1 0 ]; % final conds
% constants.N = 70;
% constants.obstacles_circles = [ 5, 5, 3];
% % constants.obstacles = [ 5 5 ] + sqrt(2)/2*[-1 1; 1 1; 1 -1; -1 -1];
% constants.statebounds = [
%     -Inf, -Inf, -2*pi, 0, -1; % inferior bounds of each state variable
%     Inf, Inf, 2*pi, 5, 1; % superior bounds of each state variable
% ];
% constants.uselogbar = true;


%%%%%%%%%%%%%%%% Example 5 %%%%%%%%%%%%%%%%
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

%%%%%%%%%%%%%%%% Example no constraints  %%%%%%%%%%%%%%%%

constants.T = 60;
constants.xi = [ 0 0 0 1 0];
constants.xf = [ 30 30 pi/2 1 0];
constants.obstacles_circles = [18,10,8];
% constants.obstacles = [18,10] + [1 1; 1 -1; -1 -1; -1 1]*8/sqrt(2);
constants.N = 20;
constants.uselogbar = true;
constants.statebounds = [
    -Inf, -Inf, -2*pi, 0, -pi/4; % inferior bounds of each state variable
    Inf, Inf, 2*pi, 1.1, pi/4; % superior bounds of each state variable
];

%%%%%%%%%%%%%%%% Diagonal Example  %%%%%%%%%%%%%%%%
% constants.T = 80;
% constants.xi = [ 0 0 pi/4 1 0];
% constants.xf = [ 30 30 pi/4 1 0];
% % constants.obstacles_circles = [15,10,8];
% % constants.obstacles = [18,10] + [1 1; 1 -1; -1 -1; -1 1]*8/sqrt(2);
% constants.N = 10;
% % constants.uselogbar = true;
% constants.statebounds = [
%     -Inf, -Inf, -2*pi, 0, -pi/4; % inferior bounds of each state variable
%     Inf, Inf, 2*pi, 1.1, pi/4; % superior bounds of each state variable
% ];
% constants.plotboatsize=2.5;

%%%%%%%%%%%%%%%% Functions %%%%%%%%%%%%%%%%
constants.costfun_single = @costfun_single;
constants.dynamics = @dynamicsDubin;
constants.init_guess = @init_guess;
constants.recoverxy = @recoverplot;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[xOut,JOut] = run_problem(constants); % runs once for a desired N
% [xOut,JOut, X_history, J_history, Times_history] = run_problem_progressive_n(constants); % runs with increasing N

%%

disp(['The final cost is ', num2str(JOut)])
%constraint_evaluation(xOut,constants);
% figure, 
% axis equal, hold on, axis ij,camroll(90)
% figure
plot_xy(xOut,constants);
legend('Bernstein x,y Plot');%, 'Vehicle position in uniform \Delta t')
xlabel('y'), ylabel('x')
% txt = string(['Runtime: ','']);
% annotation('textbox','String',txt);
plotedit('on');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t,xy, tenpoints] = recoverplot(X,constants)

    v = @(t) BernsteinEval(X(:,4),constants.T,t);
    w = @(t) BernsteinEval(X(:,5),constants.T,t);

    [t,xy] = ode45(@(t,xy)odefunc(t,xy,v,w), [0 constants.T], X(1,1:3));

    tenpoints =  BernsteinEval(X,constants.T,linspace(0,constants.T,10));

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
    u = X(:,4);
    w = X(:,5);

    ceq = [
        (DiffMat*xp) - u.*cos(psi)
        (DiffMat*yp) - u.*sin(psi)
        DiffMat*psi - w
    ];
    
    c = [];

end


function J = costfun_single(X,constants)
    % the matrix X's columns have, control points for x,y,psi,v,w variables
    % for 1 vehicle, respectfully

    v = X(:,4);
    J = sum(v.^2)/constants.N;

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
