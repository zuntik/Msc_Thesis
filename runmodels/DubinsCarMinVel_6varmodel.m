clear; close all;

addpath('..\Bernstein');
addpath('..\BeBOT_lib');

% Settings

CONSTANTS.T = 10; % time interval
CONSTANTS.xi = [5 3 cos(0) sin(0) 1 0]; % x y c s v w
CONSTANTS.xf = [0 0 cos(0) sin(0) 1 0]; % x y c s v w
CONSTANTS.N = 20;
CONSTANTS.DiffMat = BernsteinDerivElevMat(CONSTANTS.N,CONSTANTS.T);
CONSTANTS.numvars = length(CONSTANTS.xi);
CONSTANTS.obstacle = [];

xIn = init_guess(CONSTANTS);

options = optimoptions(@fmincon,'Algorithm','sqp','Display','Iter','MaxFunctionEvaluations',300000,'StepTolerance',eps,'MaxIterations',Inf);

xOut = fmincon(@(x)costfun(x,CONSTANTS),xIn,[],[],[],[],[],[],@(x)nonlcon(x,CONSTANTS),options);
xOut = [CONSTANTS.xi; reshape(xOut,[],CONSTANTS.numvars); CONSTANTS.xf];
%%
figure
BernsteinPlot(xOut(:,1:2),CONSTANTS.T);
[~,xy] = recoverplot(xOut,CONSTANTS.T);
plot(xy(:,1),xy(:,2));

figure
c = xOut(:,3);
s = xOut(:,4);
plot(c.^2+s.^2);


function [t,xy] = recoverplot(X,T)

    v = @(t) BernsteinEval(X(:,5),T,t);
    w = @(t) BernsteinEval(X(:,6),T,t);
    
    [t,xy] = ode45(@(t,xy)odefunc(t,xy,v,w), [0 T], X(1,1:4));

    function dydt = odefunc(t,y,v,w)

        dydt = zeros(4,1);

        dydt(1) =  v(t)*y(3); % dx
        dydt(2) =  v(t)*y(4); % dy
        dydt(3) = -w(t)*y(4); % dc
        dydt(4) =  w(t)*y(3); % ds

    end

end

function [c,ceq] = nonlcon(X,CONSTANTS)

    DiffMat = CONSTANTS.DiffMat;
    
    X = [CONSTANTS.xi; reshape(X,[],CONSTANTS.numvars); CONSTANTS.xf];
    xp = X(:,1);
    yp = X(:,2);
    c = X(:,3);
    s = X(:,4);
    v = X(:,5);
    w = X(:,6);

    ceq = [
        DiffMat*xp - v.*c
        DiffMat*yp - v.*s
        DiffMat*c + s.*w
        DiffMat*s - c.*w
        c.^2 + s.^2 - 1
        %v.^2 - (DiffMat*xp).^2 - (DiffMat*xp).^2
    ];

    c=-v+0.2;

end

function J = costfun(X,CONSTANTS)

    X = reshape(X,[],CONSTANTS.numvars);
    X = [CONSTANTS.xi; X; CONSTANTS.xf];
    v = X(:,5);
    w = X(:,6);
    a = CONSTANTS.DiffMat*v;
    J = sum(a.^2)+2*sum(w.^2);
    
end

function xinit = init_guess(CONSTANTS)

    N = CONSTANTS.N; 

    x = rand(N-1,1);
    y = rand(N-1,1);
    c = rand(N-1,1);
    s = sqrt(1-c.^2);
    v = rand(N-1,1);
    w = rand(N-1,1);

    xinit = [x;y;c;s;v;w];
    
end
