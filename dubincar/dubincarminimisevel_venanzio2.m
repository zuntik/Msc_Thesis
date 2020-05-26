clear; close all;

global CONSTANTS

addpath('..\Bernstein');
addpath('..\BeBOT_lib');

%% Settings
CONSTANTS.T = 10; % time interval

%% Boundary Conditions
psi_0 = 0;
psi_f = 0;

CONSTANTS.init_conds  = [5 3 psi_0 1 0]; % 
CONSTANTS.final_conds = [0 0 psi_f 1 0]; % 

%% Discretization
CONSTANTS.N = 20;


%% Initial guess
x_init = init_guess(CONSTANTS);




%% Run
options = optimoptions(@fmincon,'Algorithm','sqp','MaxFunctionEvaluations',300000,'StepTolerance',eps,'MaxIterations',Inf);
[x,f] = fmincon(@(x)costfun(x,CONSTANTS),x_init,[],[],[],[],[],[],@(x)nonlcon(x,CONSTANTS),options);






%% Plot
N = CONSTANTS.N;

xp = [CONSTANTS.init_conds(1) x(1:N-1)' CONSTANTS.final_conds(1)];
yp = [CONSTANTS.init_conds(2) x(N:2*N-2)' CONSTANTS.final_conds(2)];
psi = [CONSTANTS.init_conds(3) x(2*N-1:3*N-3)' CONSTANTS.final_conds(3)];
v = [CONSTANTS.init_conds(4) x(3*N-2:4*N-4)' CONSTANTS.final_conds(4)];
omega = [CONSTANTS.init_conds(5) x(4*N-3:5*N-5)' CONSTANTS.final_conds(5)];


t = [0:0.001:CONSTANTS.T];

figure, hold on
plot(BernsteinPoly(xp,t),BernsteinPoly(yp,t)); hold on

[~,xy] = recoverplot([xp;yp;psi;v;omega].',CONSTANTS.T);
plot(xy(:,1),xy(:,2));

% Run simulation model
out = sim('runsim');
plot(out.x.signals.values,out.y.signals.values)


legend('result','recovered','simulink_recovered')




%% Functions

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
    N = CONSTANTS.N; 
    T = CONSTANTS.T;

    xp = [CONSTANTS.init_conds(1) x(1:N-1)' CONSTANTS.final_conds(1)];
    yp = [CONSTANTS.init_conds(2) x(N:2*N-2)' CONSTANTS.final_conds(2)];
    psi = [CONSTANTS.init_conds(3) x(2*N-1:3*N-3)' CONSTANTS.final_conds(3)];
    v = [CONSTANTS.init_conds(4) x(3*N-2:4*N-4)' CONSTANTS.final_conds(4)];
    omega = [CONSTANTS.init_conds(5) x(4*N-3:5*N-5)' CONSTANTS.final_conds(5)];

    [~,~,Diff] = BeBOT(N,T);

    dyn1 = (xp*Diff-v.*cos(psi))';
    dyn2 = (yp*Diff-v.*sin(psi))';
    dyn3 = (psi*Diff - omega)';

    c=[-v'];
    ceq=[
        dyn1; dyn2; dyn3; ...
        ];
end







function J = costfun(x,CONSTANTS)
%COSTFUN Summary of this function goes here
    N = CONSTANTS.N; 
    T = CONSTANTS.T; 

    xp = [CONSTANTS.init_conds(1) x(1:N-1)' CONSTANTS.final_conds(1)];
    yp = [CONSTANTS.init_conds(2) x(N:2*N-2)' CONSTANTS.final_conds(2)];
    psi = [CONSTANTS.init_conds(3) x(2*N-1:3*N-3)' CONSTANTS.final_conds(3)];
    v = [CONSTANTS.init_conds(4) x(3*N-2:4*N-4)' CONSTANTS.final_conds(4)];
    omega = [CONSTANTS.init_conds(5) x(4*N-3:5*N-5)' CONSTANTS.final_conds(5)];

    J = 0.1*sum(v.^2)+20*sum(omega.^2);
end



function xinit = init_guess(CONSTANTS)
%COSTFUN Summary of this function goes here
    N = CONSTANTS.N; 
    T = CONSTANTS.T; 

    [~,~,Diff] = BeBOT(N,T);


    xp = 5*rand(N-1,1);
    yp = 5*rand(N-1,1);
    psi = rand(N-1,1);
    v = ones(N-1,1);
    omega = zeros(N-1,1);

    xinit = [xp;yp;psi;v;omega];
end

