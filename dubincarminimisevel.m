clear; close all;

addpath('Bernstein');
addpath('BeBOT_lib');

%% Settings

T = 10; % time interval

%% Boundary Conditions
ang_0 = 0;
ang_f = 0;

init_conds  = [3 5 cos(ang_0) sin(ang_0) 2 0]; % x y r1 r2 v w
final_conds = [0 0 cos(ang_f) sin(ang_f) 0 0]; % inf <=> don't care

%% run

N = [ 4 5 8 10 13 20 40 50 100 ];
N = 5:12;
% N = logspace(log10(5),log10(100),8);
% N = [7 10 14];
N = 14;

sols = cell(1,length(N));

%shapes = [  1 1.5 ; 1 2.5 ; 2 2.5 ; 2 1.5 ];
shapes = [];
nonlcon = @(X,N,T) nonlcon_obstacle(X,N,T,shapes);
%nonlcon = @basicnonlcon;

for i = 1:length(N)
    [xIn,xOut,dur,Jout,constraints] = run_problem(init_conds,final_conds,N(i),T,nonlcon,@costFunc);
    sols{i} = struct('N',N(i),'dur',dur,'xIn',xIn,'xOut',xOut,'Jout',Jout,'constraints',constraints);
end

%% post process

if length(N) ~= 1
    figure(1), sgtitle('x y plane')
    figure(2), sgtitle('compare vels')
    figure(3), sgtitle('compare angles')
    for i = 1:length(N)
        s = sols{i};
        plotstuff(s.xIn,s.xOut,T,shapes,2,4,i)
    end
    evolution(sols,N);
else
    figure(1), title('x y plane')
    figure(2), title('compare vels')
    figure(3), title('compare angles')
    plotstuff(xIn,xOut,T,shapes)
    testconstraints(xIn, xOut, N, T, constraints)
end

%% Plot stuff
function plotstuff(xIn,xOut, T,shapes,m,n,i)

    N = size(xIn,1)-1;

    dxOut = BernsteinDeriv(xOut,T);
    dxOut_squared = BernsteinPow(dxOut,2);

%     figure(4), grid on, axis equal
%     title('x y plane initial guess')
%     BernsteinPlot(xIn(:,1:2),T,'PlotControlPoints',false);

    figure(1)
    if nargin > 4
        subplot(m,n,i)
        title(['N = ',num2str(N)])
    end
    grid on, axis equal
    xlabel('x (m)')
    ylabel('y (m)')
    BernsteinPlot(xOut(:,1:2),T,'PlotControlPoints',false);
    [~,xy] = recoverplot(xOut,T);
    plot(xy(:,1),xy(:,2));
    if ~isempty(shapes)
        for i = 1:size(shapes,3)
            plot(polyshape(shapes(:,:,i)))
        end
    end
    legend('result','recovered')
    
    figure(2)
    if nargin > 4
        subplot(m,n,i)
        title(['N = ',num2str(N)])
    end
    grid on, hold on
    xlabel('t (s)');
    ylabel('vel (m/s)');
    BernsteinPlot(xOut(:,5),T,'PlotControlPoints',false);
    vel = @(times)arrayfun(@(t)sqrt(sum(BernsteinEval(dxOut_squared(:,1:2),T,t),2)),times);
    fplot(vel,[0 T]);
    legend('v','sqrt(dx^2+dy^2)')

    figure(3)
    if nargin > 4
        subplot(m,n,i)
        title(['N = ',num2str(N)])
    end
    grid on, hold on
    xlabel('t (s)')
    ylabel('ang (rad)')
    ang  = @(times)arrayfun(@(t)atan2(BernsteinEval(xOut(:,4),T,t),BernsteinEval(xOut(:,3),T,t)),times);
    ang2 = @(times)arrayfun(@(t)atan2(BernsteinEval(dxOut(:,2),T,t),BernsteinEval(dxOut(:,1),T,t)),times);
    fplot(ang,[0 T]);
    fplot(@(t)ang2(t),[0 T]);
    legend('atan2(c,s)','atan2(dy,dx)');

end

function testconstraints(xIn, xOut, N, T, constraints)

    dxOut = BernsteinDeriv(xOut,T);
    dxOut_squared = BernsteinPow(dxOut,2);
    
    disp(['initial guess respects upper bounds? ',num2str(all(reshape(constraints.ub,[N+1 6])-xIn>=0,'all'))])
    disp(['initial guess respects lower bounds? ',num2str(all(xIn-reshape(constraints.lb,[N+1 6]) >=0,'all'))])

    disp(['norm of the linear constraints error: ', num2str(norm(constraints.Aeq*xOut(:)-constraints.beq))]);
    [~,ceq]=constraints.nonlcon(xOut);
    disp(['norm of the nonlinear constraints error: ', num2str(norm(ceq))]);

    vel = @(times)arrayfun(@(t)sqrt(sum(BernsteinEval(dxOut_squared(:,1:2),T,t),2)),times);
    disp(['first true vel: ',num2str(vel(0))]);
    disp(['last true vel: ',num2str(vel(T))]);

    ang = @(times)arrayfun(@(t)atan2(BernsteinEval(xOut(:,4),T,t),BernsteinEval(xOut(:,3),T,t)),times);
    disp(['first ang: ',num2str(ang(0))]);
    disp(['last ang: ',num2str(ang(T))]);

    disp(['sol respect upper bounds? ',num2str(all(reshape(constraints.ub,[N+1 6])-xOut>=0,'all'))]);
    disp(['sol respect lower bounds? ',num2str(all(xOut-reshape(constraints.lb,[N+1 6]) >=0,'all'))]);

end

function evolution(sols,N)

    figure(4)
    %set(gcf,'WindowState','fullscreen')

    subplot(2,2,1)
    plot(N,arrayfun(@(i)sols{i}.Jout,1:length(N)))
    title('cost')
    xlabel('N (order)'), ylabel('Jout')

    subplot(2,2,2)
    plot(N,arrayfun(@(i)sols{i}.dur,1:length(N)))
    title('duration')
    xlabel('N (order)'), ylabel('runtime (s)')

    subplot(2,2,3)
    plot(N,arrayfun(@(i)norm(sols{i}.constraints.Aeq*sols{i}.xOut(:)-sols{i}.constraints.beq),1:length(N)));
    title('norm linear error')
    xlabel('N (order)'), ylabel('linear error')

    subplot(2,2,4)
    nonlcon_norms = zeros(1,length(N)); 
    for i = 1:length(N)
        [~,ceq] = sols{i}.constraints.nonlcon(sols{i}.xOut(:));
        nonlcon_norms(i) = norm(ceq);
    end
    plot(N,nonlcon_norms);
    title('norm nonlinear error')
    xlabel('N (order)'), ylabel('nonlinear error')
end


%% run problem
function [xIn,xOut,dur,Jout,constraints] = run_problem(xi,xf,N,T,dynamics,costFunc)
    
%     Nv = size(xi,1);
%     if Nv  ~= 1
%         disp('no support for multiple vehicles yet');
%     end

    [Aeq,beq,lb,ub] = LinConstr(N,T,xi,xf);

    nonlcon=@(x) dynamics(x,N,T);
    
    xIn = InitialGuess(N, T, xi, xf);
    
    % options = optimoptions(@fmincon,'Algorithm','sqp',...
    %                        'MaxFunctionEvaluations',40000,...
    %                        'ConstraintTolerance',1e-6,...
    %                        'StepTolerance',1e-6,...
    %                        'Display','iter');
    % %                        'OptimalityTolerance',1e-8,...
    % % %                       'FunctionTolerance',1e-8,..
    % %                       'ConstraintTolerance',1e-3);

    options = optimoptions(@fmincon,'Algorithm','sqp');

    tic;
        [xOut,Jout,~,~] = fmincon(@(x)costFunc(x,N,T),xIn,[],[],Aeq,beq,lb,ub,nonlcon,options);
    dur=toc;
    
    constraints=struct('Aeq',Aeq,'beq',beq,'lb',lb,'ub',ub,'nonlcon',nonlcon);
    
end

%% Cost function
function J = costFunc(X,N,T)

    [~, ~, ~, ~, ~, ~, ~, ~, v, w]  = getControlPoints(X,N,T);

    b_sq_int = @(p) BernsteinIntegr(BernsteinPow(p,2),T);

    cw = 1;
    cv = 1;
    
    J = cw * b_sq_int(w) + cv* b_sq_int(v);

end

function [t,xy] = recoverplot(X,T)

    v = @(t) BernsteinEval(X(:,5),T,t);
    w = @(t) BernsteinEval(X(:,6),T,t);
    
    [t,xy] = ode45(@(t,xy)odefunc(t,xy,v,w), [0 T], [X(1,1) X(1,2) X(1,3) X(1,4)]);

    function dydt = odefunc(t,y,v,w)

        dydt = zeros(4,1);

        dydt(1) =  v(t)*y(3);
        dydt(2) =  v(t)*y(4);
        dydt(3) = -w(t)*y(4);
        dydt(4) =  w(t)*y(3);

    end

end

%% Linear contraints
% The linear contraints set the boundary conditions
function [Aeq, beq,lb,ub] = LinConstr(N,T,xi,xf)

    % initial conds
    Aeq = zeros(6,(N+1)*6);
    for i = 1:6
        Aeq(i,1+(i-1)*(N+1)) = 1;
    end
    beq = xi(:);

    % final conds
    xf=xf(xf~=Inf);
    Aeq=[Aeq;zeros(length(xf),(N+1)*6)];
    for i = 1:length(xf)
        Aeq(i+6, i*(N+1)) = 1;
    end
    beq=[beq;xf(:)];

    % optional part: help with derivatives
    if true
        Aeq(end+1,1:2) = [ -N/T, N/T ];
        Aeq(end+1, (N+1)+1:(N+1)+2) = [ -N/T, N/T ];
        beq = [beq ; xi(5).*xi(3:4).'];

        Aeq(end+1, (N+1)-1:(N+1)) = [ -N/T, N/T ];
        Aeq(end+1, (N+1)*2-1:(N+1)*2) = [ -N/T, N/T ];
        beq = [beq ; xf(5).*xf(3:4).'];
    end
    
    if xi(5) == 0 && true
        disp('using cool constraint for initial condition')
        Aeq(end+1,3) = xi(4);
        Aeq(end,(N+1)+3) = xi(3);
        beq = [beq; xi(4)*xi(1)-xi(3)*xi(2)];
    end
    if xf(5) == 0 && true
        disp('using cool constraint for final condition')
        Aeq(end+1,(N+1)-2) = xf(4);
        Aeq(end,2*(N+1)-2) = -xf(3);
        beq = [beq; xf(4)*xf(1)-xf(3)*xf(2)];
    end
    
    % variable bounds

    % states
    x_ub = ones(4*(N+1),1)*Inf;
    x_lb = -x_ub;

    % control 

    omg_max = 10*180/pi;
    vmax = norm(xf(1:2)-xi(1:2))/T*1.5;
    vmax = max([xi(5), xf(5), vmax]);
    u_ub = ones((N+1),2).*[vmax, omg_max];
    u_lb = ones((N+1),2).*[0, -omg_max]; 

    lb = [x_lb(:); u_lb(:)];
    ub = [x_ub(:); u_ub(:)];

end

%% Nonlinear contraints
% The nonlinear contraints set the dynamic conditions
% this function is named "basic" because it doesn't take obstacles into
% account
function [c, ceq] = basicnonlcon(X,N,T)

    [~, ~, c, s, dx, dy, dc, ds, v, w] = getTrajectories(X,N,T);

    ceq = [
        dx - v.*c
        dy - v.*s
        dc + s.*w
        ds - c.*w
        c.^2 + s.^2 - 1
        %v.^2 - dx.^2 - dy.^2
    ];

    c = [];

end

function [c, ceq] = nonlcon_obstacle(X,N,T,shapes)
    [~, ceq] = basicnonlcon(X,N,T);

    if ~isempty(shapes)
        c = zeros(size(shapes,3),1);
        for i = 1:size(shapes,3)
            c = 0.01-MinDistBernstein2Polygon(X(:,[1:2]).', shapes(:,:,i).');
        end
    else
        c = [];
    end

end


%% Initial Guess
function [X] = InitialGuess(N, T, init_conds, final_conds)
    X = zeros(N+1,6);

    final_conds(final_conds==Inf) = 0;

    vi = init_conds(5)  .* init_conds(3:4);
    vf = final_conds(5) .* final_conds(3:4);

    X(1,1:6) = init_conds;
    X(end,1:6) = final_conds;

    min_dist = norm(final_conds(1:2)-init_conds(1:2));

    % deal with the second and penultimate control points
    X(2,1:2) = (T/N).*vi + init_conds(1:2);
    X(end-1,1:2) = -(T/N).*vf + final_conds(1:2);

    % deal with third and antepenultimate
    X(3,1:2) =     X(2,1:2)     + (min_dist/N).* init_conds(3:4);
    X(end-2,1:2) = X(end-1,1:2) - (min_dist/N).*final_conds(3:4);
    
    % place remaning control points in a straight line
%     X(3:end-2,1) = linspace(X(3,1),X(end-2,1),N-3);
%     X(3:end-2,2) = linspace(X(3,2),X(end-2,2),N-3);
    
    %%%%%
    % this code is to avoid the obstacle
%     mid = floor(N/2);
%     X(mid,1:2) = [0 3];
%     X(3:mid,1) = linspace(X(3,1),X(mid,1),size(X(3:mid,1),1));
%     X(3:mid,2) = linspace(X(3,2),X(mid,2),size(X(3:mid,1),1));
%     X(mid:end-2,1) = linspace(X(mid,1),X(end-2,1),size(X(mid:end-2,1),1));
%     X(mid:end-2,2) = linspace(X(mid,2),X(end-2,2),size(X(mid:end-2,1),1));
    %%%%%
    
    vec = X(end-2,1:2) - X(3,1:2);
    
    ang_mid = atan2(vec(2),vec(1));

    for it = 3:(size(X,1)-3)
        X(it,3:4) = [ cos(ang_mid) sin(ang_mid) ];
    end
    X(2,3:4) = init_conds(3:4);
    vec_fin = X(end,1:2)-X(end-2,1:2);
    for it = (-2:0)+size(X,1)
        ang_fin = atan2(vec_fin(2),vec_fin(1));
        X(it,3:4) = [ cos(ang_fin) sin(ang_fin) ];
    end
    
end

%% Get Control Points of Trajectories
function [px, py, c, s, dx, dy, dc, ds, v, w] = getControlPoints(X,N,T)

    persistent old_N old_T d_mat
    if  isempty(old_N) || isempty(old_T) || old_N ~= N || old_T ~= T
        old_N = N;
        old_T = T;
        d_mat = BernsteinDerivElevMat(N,T);
    end

    X = reshape(X(:), [N+1,6]); % this step could be redundant

    % get control points
    px = X(:,1);
    py = X(:,2);
    c = X(:,3);
    s = X(:,4);
    v = X(:,5);
    w = X(:,6);

    % compute  derivatives 
    dX  = d_mat*X;
    dx  = dX(:,1);
    dy  = dX(:,2);
    dc = dX(:,3);
    ds = dX(:,4);

end

%% Get Trajectories
function [px, py, c, s, dx, dy, dc, ds, v, w] = getTrajectories(X,N,T)

    persistent old_N old_T d_mat eval_mat
    if  isempty(old_N) || isempty(old_T) || old_N ~= N || old_T ~= T
        old_N = N;
        old_T = T;
        d_mat = BernsteinDerivElevMat(N,T);
        eval_mat = BernsteinCtrlPntsEval(N);
    end

    X = reshape(X(:), [N+1,6]); % this step could be redundant

    dX = d_mat*X;
    X  = eval_mat*X;
    dX = eval_mat*dX;

    % get control points
    px = X(:,1);
    py = X(:,2);
    c = X(:,3);
    s = X(:,4);
    v = X(:,5);
    w = X(:,6);

    % compute  derivatives 
    dx  = dX(:,1);
    dy  = dX(:,2);
    dc = dX(:,3);
    ds = dX(:,4);

end