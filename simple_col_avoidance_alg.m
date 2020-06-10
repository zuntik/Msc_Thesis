clear all; close all;

addpath('.\Bernstein');
addpath('.\BeBOT_lib');

% Settings
T = 10;
N = 13;

% Boundary Conditions

xi  = [3 5 1 0]; % x y v w
xf = [0 0 1 0]; % inf <=> don't care


shapes = [  1 1.5 ; 1 2.5 ; 2 2.5 ; 2 1.5 ];
%nonlcon = @(X,N,T) nonlcon_obstacle(X,N,T,shapes);
nonlcon = @basicnonlcon;


xOut = run_problem(xi,xf,N,T,nonlcon,@costFunc);



% run problem
function xOut = run_problem(xi,xf,N,T,dynamics,costFunc)

    nonlcon=@(x) dynamics(x,xi,xf,N,T);
    
    xIn = InitialGuess(N);

    %options = optimoptions(@fmincon,'Algorithm','sqp','Display','iter','MaxFunctionEvaluations',1000*(N-1)*6);
    options = optimoptions(@fmincon,'Display','Iter','Algorithm','sqp','MaxFunctionEvaluations',300000,'StepTolerance',eps,'MaxIterations',Inf);
    
    xOut = fmincon(@(x)costFunc(x,xi,xf,N,T),xIn,[],[],[],[],[],[],nonlcon,options);
    
    xOut = [xi;reshape(xOut(:),[],4);xf];
        
end

% Cost function
function J = costFunc(X,xi,xf,N,T)
    [~,~,v,w]  = getControlPoints(X,xi,xf);
    a = BernsteinDeriv(v,T);
    ca = 1;
    cw = 2;
    J = ca*sum(a.^2)+cw*sum(w.^2);
end


% Nonlinear contraints
% The nonlinear contraints set the dynamic conditions
% this function is named "basic" because it doesn't take obstacles into
% account
function [c, ceq] = basicnonlcon(X,xi,xf,N,T)

    persistent old_N old_T degr_elev_mat degr_elev_mat2
    if  isempty(old_N) || isempty(old_T) || old_N ~= N || old_T ~= T
        old_N = N;
        old_T = T;
        degr_elev_mat = BernsteinDegrElevMat(2*N-3,3*N-2);
        degr_elev_mat2 = BernsteinDegrElevMat(2*N-2,2*N);
    end


    [x,y,v,w] = getControlPoints(X,xi,xf);
    
    dx = BernsteinDeriv(x,T);
    dy = BernsteinDeriv(y,T);
    ddx = BernsteinDeriv(dx,T);
    ddy = BernsteinDeriv(dy,T);
    dx_2 = BernsteinPow(dx,2);
    dy_2 = BernsteinPow(dy,2);
    
    left = BernsteinMul(dx_2 + dy_2,w);
    right = BernsteinMul(dx,ddy)-BernsteinMul(ddx,dy);
    
    left2 = BernsteinPow(v,2);
    right2 = dx_2 + dy_2;
    
    ceq = [
%         left-BernsteinDegrElev(right,3*N-2)
%         left2-BernsteinDegrElev(right2,2*N);
        left-degr_elev_mat*right
        left2-degr_elev_mat2*right2
        v(1) - xi(3)
        w(1) - xi(4)
        v(end) - xf(3)
        w(end) - xf(4)
    ];
    
    c = -v;

end

function [c, ceq] = nonlcon_obstacle(X,xi,xf,N,T,shapes)
    [~, ceq] = basicnonlcon(X,xi,xf,N,T);

    if ~isempty(shapes)
        c = zeros(size(shapes,3),1);
        for i = 1:size(shapes,3)
            c = 0.01-MinDistBernstein2Polygon(X(:,[1:2]).', shapes(:,:,i).');
        end
    else
        c = [];
    end

end

% Initial Guess
function [X] = InitialGuess(N)

    x = rand(N-1,1);
    y = rand(N-1,1);
    v = rand(N-1,1);
    w = rand(N-1,1);
    X = [ x;y;v;w ];

end

function [x,y,v,w] = getControlPoints(X,xi,xf)
    X = reshape(X(:), [], 4);
    X = [ xi; X; xf ];
    x = X(:,1);
    y = X(:,2);
    v = X(:,3);
    w = X(:,4);    
end
