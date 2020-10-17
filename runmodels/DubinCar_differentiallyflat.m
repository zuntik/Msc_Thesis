clear all; %#ok<CLALL>

addpath('..\Bernstein');
addpath('..\BeBOT_lib');
addpath('..\TrajecOptimLib');

constants.T = 10;
constants.xi = [ 0 0]; % 0 1 0];
constants.xf = [ 5 5]; % pi/2 1 0];
constants.N = 30;
constants.obstacles_circles = [5,0,3];

% extra stuff
constants.psii=0;
constants.psif=pi/2;
constants.vi=1;
constants.vf=1;
constants.wi=0;
constants.wf=0;
constants.times = 0:constants.T/1000:constants.T;
constants.evalmat_anum=BernsteinEvalMat(constants.N*2-3,constants.T,constants.times);
constants.evalmat_adenumsquare=BernsteinEvalMat(constants.N*2-2,constants.T,constants.times);

% functions
constants.costfun_single = @costfun_single;
constants.dynamics = @dynamicsDubin;
% constants.init_guess = @init_guess;
constants.recoverxy = @recoverplot;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[xOut,JOut] = run_problem(constants);

disp(['The final cost is ', num2str(JOut)])

BernsteinPlot(xOut, constants.T);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [c,ceq] = dynamicsDubin(X,constants)

    xp = X(:,1);
    yp = X(:,2);
    N = constants.N;
    T = constants.T;

    ceq = [
        % v
        constants.vi-sqrt(((xp(2)-xp(1))*N/T)^2+((yp(2)-yp(1))*N/T)^2)
        constants.vf-sqrt(((xp(N+1)-xp(N))*N/T)^2+((yp(N+1)-yp(N))*N/T)^2)
        % psi
        constants.psii-atan2(yp(2)-yp(1),xp(2)-xp(1))
        constants.psif-atan2(yp(N+1)-yp(N),xp(N+1)-xp(N))
    ];

    c = [];
    
end

function J = costfun_single(X, constants)

    xp = X(:,1);
    yp = X(:,2);
    T = constants.T;
    
    dx = BernsteinDeriv(xp, T);
    ddx = BernsteinDeriv(dx, T);
    dy = BernsteinDeriv(yp, T);
    ddy = BernsteinDeriv(dy, T);
    
    a_num = BernsteinMul(dx,ddx)+BernsteinMul(dy,ddy);
    a_denum_square = eps+BernsteinPow(dx,2)+BernsteinPow(dy,2);
    
    a_num_vals = constants.evalmat_anum*a_num;
    a_denum_vals = sqrt(constants.evalmat_adenumsquare*a_denum_square);
    
    r_num_vals = constants.evalmat_anum*(BernsteinMul(dx,ddy)-BernsteinMul(ddx,dy));
    r_denum_vals = constants.evalmat_adenumsquare*(eps+BernsteinPow(dx,2)+BernsteinPow(dy,2));
    
    a_vals = a_num_vals./a_denum_vals;
    r_vals = r_num_vals./r_denum_vals;
    J = trapz(constants.times,a_vals.^2)+2*trapz(constants.times,r_vals.^2);
    
end

function J = costfun_single_slow(X, constants) %#ok<DEFNU>

    xp = X(:,1);
    yp = X(:,2);
    T = constants.T;
    
    dx = BernsteinDeriv(xp, T);
    ddx = BernsteinDeriv(dx, T);
    dy = BernsteinDeriv(yp, T);
    ddy = BernsteinDeriv(dy, T);
    a_num=@(t)BernsteinEval(BernsteinMul(dx,ddx)+BernsteinMul(dy,ddy),T,t);
    a_denum=@(t)sqrt(BernsteinEval(eps+BernsteinPow(dx,2)+BernsteinPow(dy,2),T,t));
    r_num=@(t)BernsteinEval(BernsteinMul(dx,ddy)-BernsteinMul(ddx,dy),T,t);
    r_denum=@(t)BernsteinEval(eps+BernsteinPow(dx,2)+BernsteinPow(dy,2),T,t);
    
    a = @(t) a_num(t)/a_denum(t);
    asquared = @(t) a(t)^2;
    
    r = @(t) r_num(t)/r_denum(t);
    rsquared = @(t) r(t)^2;
    
    asquared = @(t) arrayfun(asquared,t);
    rsquared = @(t) arrayfun(rsquared,t);
    
    J = integral(asquared, 0, T) + integral(rsquared, 0, T);
    
end

function xinit = init_guess(constants) %#ok<DEFNU>

    N = constants.N; 

    xinit = zeros((constants.N-1)*constants.numvars,constants.Nv);
    for i = 1:constants.Nv
        x = linspace(constants.xi(i,1),constants.xf(i,1),N-1).';
        y = linspace(constants.xi(i,2),constants.xf(i,2),N-1).';
        xinit(:,i) = [x;y];
    end
    %xinit = xinit(:);
    
end
