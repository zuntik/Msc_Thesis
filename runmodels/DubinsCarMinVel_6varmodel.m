clear; close all;

addpath('..\Bernstein');
addpath('..\BeBOT_lib');
addpath('..\TrajecOptimLib\');

% Settings

constants.T = 100; % time interval
constants.xi = [0 0 cos(pi/4) sin(pi/4) .1 0]; % x y c s v w
constants.xf = [10 10 cos(pi/4) sin(pi/4) .1 0]; % x y c s v w
constants.N = 20;
constants.obstacles_circles = [ 5, 5, 1];
constants.statebounds = [
    -Inf, -Inf, -1, -1, 0, -1; % inferior bounds of each state variable
    Inf, Inf, 1, 1, 5, 1; % superior bounds of each state variable
];


%%%%%%%%%%%%%%%% Functions %%%%%%%%%%%%%%%%
constants.costfun_single = @costfun_single;
constants.dynamics = @dynamicsDubin;
constants.init_guess = @init_guess;
constants.recoverxy = @recoverplot;


%%

[xOut,JOut] = run_problem(constants); % runs once for a desired N

plot_xy(xOut,constants);

figure
c = xOut(:,3);
s = xOut(:,4);
plot(c.^2+s.^2);


function [t,xy, tenpoints] = recoverplot(X,constants)

    v = @(t) BernsteinEval(X(:,5),constants.T,t);
    w = @(t) BernsteinEval(X(:,6),constants.T,t);
    
    [t,xy] = ode45(@(t,xy)odefunc(t,xy,v,w), [0 constants.T], X(1,1:4));

    tenpoints =  BernsteinEval(X,constants.T,linspace(0,constants.T,10));
    for i = 1:10
        tenpoints(i, 3) = atan2(tenpoints(i, 4), tenpoints(i, 3));
    end
    tenpoints = tenpoints(:, 1:3);
    
    function dydt = odefunc(t,y,v,w)

        dydt = zeros(4,1);

        dydt(1) =  v(t)*y(3); % dx
        dydt(2) =  v(t)*y(4); % dy
        dydt(3) = -w(t)*y(4); % dc
        dydt(4) =  w(t)*y(3); % ds

    end

end

function [c,ceq] = dynamicsDubin(X,constants)

    DiffMat = constants.DiffMat;
    
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

%     c=-v+0.2;
    c = [];

end

function J = costfun_single(X,constants)

    v = X(:,5);
    w = X(:,6);
    a = constants.DiffMat*v;
    J = sum(a.^2)+2*sum(w.^2);
    
end

function xinit = init_guess(constants)

    N = constants.N; 

    xinit = zeros((constants.N-1)*constants.numvars,constants.Nv);
    for i = 1:constants.Nv
        x = linspace(constants.xi(i,1),constants.xf(i,1),N-1).';
        y = linspace(constants.xi(i,2),constants.xf(i,2),N-1).';
        c = rand(N-1,1);
        s = sqrt(1-c.^2);
        v = ones(N-1,1);
        w = zeros(N-1,1);
        xinit(:,i) = [x;y;c;s;v;w];
    end
    %xinit = xinit(:);

end
