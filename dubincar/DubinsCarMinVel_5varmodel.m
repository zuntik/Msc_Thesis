clear all; close all;

addpath('..\Bernstein');
addpath('..\BeBOT_lib');

% Settings

% CONSTANTS.T = 10; % time interval
CONSTANTS.T = 50; % time interval
% CONSTANTS.xi = [5 3 0 1 0]; % x y psi v w
% CONSTANTS.xf = [0 0 0 1 0]; % x y psi v w
CONSTANTS.xi  = [-10 40 0 1.1 0]; % x y psi v w
CONSTANTS.xf = [0 0 0 1.2 0]; % x y psi v w
CONSTANTS.N = 40;
CONSTANTS.DiffMat = BernsteinDerivElevMat(CONSTANTS.N,CONSTANTS.T);
CONSTANTS.numvars = length(CONSTANTS.xi);
CONSTANTS.obstacle = [];

xIn = init_guess(CONSTANTS);

options = optimoptions(@fmincon,'Display','Iter','Algorithm','sqp','MaxFunctionEvaluations',300000,'StepTolerance',eps,'MaxIterations',Inf);

xOut = fmincon(@(x)costfun(x,CONSTANTS),xIn,[],[],[],[],[],[],@(x)nonlcon(x,CONSTANTS),options);
xOut = [CONSTANTS.xi; reshape(xOut,[],CONSTANTS.numvars); CONSTANTS.xf];
%%
figure
BernsteinPlot(xOut(:,1:2),CONSTANTS.T);
[~,xy] = recoverplot(xOut,CONSTANTS.T);
plot(xy(:,1),xy(:,2));



function [t,xy] = recoverplot(X,T)

    v = @(t) BernsteinEval(X(:,4),T,t);
    w = @(t) BernsteinEval(X(:,5),T,t);
    
    [t,xy] = ode45(@(t,xy)odefunc(t,xy,v,w), [0 T], X(1,1:3));

    function dydt = odefunc(t,y,v,w)

        dydt = zeros(3,1);
        dydt(1) = v(t)*cos(y(3));%x
        dydt(2) = v(t)*sin(y(3));%y
        dydt(3) = w(t);%psi
        
    end

end

function [c,ceq] = nonlcon(X,CONSTANTS)

    DiffMat = CONSTANTS.DiffMat;
    
    X = [CONSTANTS.xi; reshape(X,[],CONSTANTS.numvars); CONSTANTS.xf];
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
    c=[-v+0.8;v-1.2 ;psi-pi;-psi-pi];

end

function J = costfun(X,CONSTANTS)

    X = [CONSTANTS.xi; reshape(X,[],CONSTANTS.numvars); CONSTANTS.xf];
    v = X(:,4);
    w = X(:,5);
    a = CONSTANTS.DiffMat*v;
    J = sum(a.^2)+2*sum(w.^2);
    
end

function xinit = init_guess(CONSTANTS)

    N = CONSTANTS.N; 

    x = rand(N-1,1);
    y = rand(N-1,1);
    psi = -ones(N-1,1);
    v = ones(N-1,1);
    w = zeros(N-1,1);

    xinit = [x;y;psi;v;w];
    
end

function xinit = bad_init_guess(CONSTANTS)

    N = CONSTANTS.N; 

    x = rand(N-1,1);
    y = rand(N-1,1);
    psi = rand(N-1,1);
    v = rand(N-1,1);
    w = rand(N-1,1);

    xinit = [x;y;psi;v;w];
    
end
