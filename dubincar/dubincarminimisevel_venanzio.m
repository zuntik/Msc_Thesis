clear; close all;

global CONSTANTS

addpath('Bernstein');
addpath('BeBOT_lib');

%% Settings
CONSTANTS.T = 10; % time interval

%% Boundary Conditions
psi_0 = 0;
psi_f = 0;

CONSTANTS.init_conds  = [5 3 psi_0 1 0]; % 
CONSTANTS.final_conds = [0 0 psi_f 1 0]; % 

%% Discretization
CONSTANTS.N = 80;


%% Initial guess
x_init = init_guess(CONSTANTS);




%% Run
options = optimoptions(@fmincon,'Algorithm','sqp','MaxFunctionEvaluations',300000,'StepTolerance',eps,'MaxIterations',Inf);
[x,f] = fmincon(@(x)costfun(x,CONSTANTS),x_init,[],[],[],[],[],[],@(x)nonlcon(x,CONSTANTS),options);


%% Plot
N = CONSTANTS.N;

x = reshape(x,[],5);
x = [CONSTANTS.init_conds; x; CONSTANTS.final_conds];
xp = x(:,1);
yp = x(:,2);
psi = x(:,3);
v = x(:,4);
omega = x(:,5);


figure 
BernsteinPlot(x(:,1:2),CONSTANTS.T,'PlotControlPoints',false), hold on

[~,xy] = recoverplot(x,CONSTANTS.T);
plot(xy(:,1),xy(:,2));
legend('result','recovered')

figure, hold on
dx = BernsteinDerivElev(xp,CONSTANTS.T);
dy = BernsteinDerivElev(yp,CONSTANTS.T);
psi_xy = @(times)arrayfun(@(t)atan2(BernsteinEval(dy,CONSTANTS.T,t),BernsteinEval(dx,CONSTANTS.T,t)),times);
psi_orig = @(times)arrayfun(@(t)BernsteinEval(psi,CONSTANTS.T,t),times)-0.1;
fplot(psi_xy, [0 CONSTANTS.T]);
fplot(psi_orig, [0 CONSTANTS.T]);


%% functions

function [t,xy] = recoverplot(X,T)

    v = @(t) BernsteinEval(X(:,4),T,t);
    w = @(t) BernsteinEval(X(:,5),T,t);
    
    [t,xy] = ode45(@(t,xy)odefunc(t,xy,v,w), [0 T], [X(1,1) X(1,2) X(1,3)]);

    function dydt = odefunc(t,y,v,w)

        dydt = zeros(3,1);
        dydt(1) = v(t)*cos(y(3));%x
        dydt(2) = v(t)*sin(y(3));%y
        dydt(3) = w(t);%psi
        
    end

end



function [c,ceq] = nonlcon(x,CONSTANTS)
    persistent old_N old_T Diff

    N = CONSTANTS.N; 
    T = CONSTANTS.T;

    if  isempty(old_N) || isempty(old_T) || old_N ~= N || old_T ~= T
        old_N = N;
        old_T = T;
        Diff = BernsteinDerivElevMat(N,T);
    end
    
    x = reshape(x,[],5);
    x = [CONSTANTS.init_conds; x; CONSTANTS.final_conds];
    xp = x(:,1);
    yp = x(:,2);
    psi = x(:,3);
    v = x(:,4);
    omega = x(:,5);

    dyn1 = (Diff*xp-v.*cos(psi))';
    dyn2 = (Diff*yp-v.*sin(psi))';
    dyn3 = (Diff*psi - omega)';

%     c=-v';
%     ceq=[
%         dyn1; dyn2; dyn3; ...
%     ];

    c=[-v'+0.2;psi'-pi;-psi'-pi];
    ceq=[
        dyn1; dyn2; dyn3; ...
        ];
end


function J = costfun(x,CONSTANTS)
%COSTFUN Summary of this function goes here
    N = CONSTANTS.N; 
    T = CONSTANTS.T; 

    x = reshape(x,[],5);
    x = [CONSTANTS.init_conds; x; CONSTANTS.final_conds];
    v = x(:,4);
    omega = x(:,5);
    acc = BernsteinDerivElev(v,T);
    
    J = sum(acc.^2)+2*sum(omega.^2);
end


function xinit = init_guess(CONSTANTS)
%COSTFUN Summary of this function goes here
    N = CONSTANTS.N; 

    xp = rand(N-1,1);
    yp = rand(N-1,1);
    psi = -ones(N-1,1);
    v = ones(N-1,1);
    omega = zeros(N-1,1);

    xinit = [xp;yp;psi;v;omega];
end


