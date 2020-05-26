clear; close all;

addpath('Bernstein');
addpath('BeBOT_lib');

% Settings
T = 10;
N = 50;

% Boundary Conditions
ang_0 = 0;
ang_f = 0;

init_conds  = [3 5 cos(ang_0) sin(ang_0) 1 0]; % x y c s v w
final_conds = [0 0 cos(ang_f) sin(ang_f) 1 0]; % inf <=> don't care

% run
sols = cell(1,length(N));

%shapes = [  1 1.5 ; 1 2.5 ; 2 2.5 ; 2 1.5 ];
%nonlcon = @(X,N,T) nonlcon_obstacle(X,N,T,shapes);
nonlcon = @basicnonlcon;

for i = 1:length(N)
    [xIn,xOut,dur,Jout,constraints] = run_problem(init_conds,final_conds,N(i),T,nonlcon,@costFunc);
    sols{i} = struct('N',N(i),'dur',dur,'xIn',xIn,'xOut',xOut,'Jout',Jout,'constraints',constraints);
end

% post process

if length(N) ~= 1
    figure(1), sgtitle('x y plane')
    figure(2), sgtitle('compare vels')
    figure(3), sgtitle('compare angles')
    for i = 1:length(N)
        s = sols{i};
        plotstuff(s.xOut,N, T,shapes,2,4,i)
    end
    evolution(sols,N);
else
    figure(1), title('x y plane')
    figure(2), title('compare vels')
    figure(3), title('compare angles')
    plotstuff(xOut, N, T,shapes)
    testconstraints(xOut, T, constraints)
end

% Functions

% Plot stuff
function plotstuff(xOut, N, T,shapes,m,n,i)

    dxOut = BernsteinDeriv(xOut,T);
    dxOut_squared = BernsteinPow(dxOut,2);

    figure(1)
    if nargin > 6
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
    if nargin > 6
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
    if nargin > 6
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

function testconstraints(xOut, T, constraints)

    dxOut = BernsteinDeriv(xOut,T);
    dxOut_squared = BernsteinPow(dxOut,2);
    
    vel = @(times)arrayfun(@(t)sqrt(sum(BernsteinEval(dxOut_squared(:,1:2),T,t),2)),times);
    disp(['first true vel: ',num2str(vel(0))]);
    disp(['last true vel: ',num2str(vel(T))]);

    ang = @(times)arrayfun(@(t)atan2(BernsteinEval(xOut(:,4),T,t),BernsteinEval(xOut(:,3),T,t)),times);
    disp(['first ang: ',num2str(ang(0))]);
    disp(['last ang: ',num2str(ang(T))]);

    disp(['sol respect upper bounds? ',num2str(all(reshape(constraints.ub,[],6)-xOut(2:end-1,:)>=0,'all'))]);
    disp(['sol respect lower bounds? ',num2str(all(xOut(2:end-1,:)-reshape(constraints.lb,[],6) >=0,'all'))]);

    [~,ceq] = constraints.nonlcon(xOut(2:end-1,:));
    disp(['nonlcon error', num2str(norm(ceq))]);
    
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

% run problem
function [xIn,xOut,dur,Jout,constraints] = run_problem(xi,xf,N,T,dynamics,costFunc)

    [lb,ub] = LinConstr(xi,xf,N,T);

    nonlcon=@(x) dynamics(x,xi,xf,N,T);
    
    xIn = InitialGuess(N,T,xi,xf);

    options = optimoptions(@fmincon,'Algorithm','sqp','Display','iter','MaxFunctionEvaluations',1000*(N-1)*6);

    tic;
        [xOut,Jout,~,~] = fmincon(@(x)costFunc(x,xi,xf,N,T),xIn,[],[],[],[],[],[],nonlcon,options);
    dur=toc;
    
    xOut = [xi;reshape(xOut(:),[],6);xf];
    
    constraints=struct('lb',lb,'ub',ub,'nonlcon',nonlcon);
    
end

% Cost function
function J = costFunc(X,init_conds,final_conds,N,T)

    [~, ~, ~, ~, ~, ~, ~, ~, v, w]  = getControlPoints(X,init_conds,final_conds,N,T);

    %b_sq_int = @(p) BernsteinIntegr(BernsteinPow(p,2),T);

    cv = 0.1;
    cw = 20;
    
    %J = cw * b_sq_int(w) + cv* b_sq_int(v);
    J = cv*sum(v.^2)+cw*sum(w.^2);
    

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

% Linear contraints
% The linear contraints set the boundary conditions
function [lb,ub] = LinConstr(xi,xf,N,T)

    % variable bounds

    % states
    x_ub = ones((N-1),4)*Inf;
    x_lb = -x_ub;

    % control 

    omg_max = 50*pi/180;
    vmax = norm(xf(1:2)-xi(1:2))/T*1.5;
    vmax = max([xi(5), xf(5), vmax]);
    u_ub = ones((N-1),2).*[vmax, omg_max];
    u_lb = ones((N-1),2).*[0, -omg_max]; 

    lb = [x_lb(:); u_lb(:)];
    ub = [x_ub(:); u_ub(:)];

end

% Nonlinear contraints
% The nonlinear contraints set the dynamic conditions
% this function is named "basic" because it doesn't take obstacles into
% account
function [c, ceq] = basicnonlcon(X,init_conds,final_conds,N,T)

    [~, ~, c, s, dx, dy, dc, ds, v, w] = getControlPoints(X,init_conds,final_conds,N,T);
    
    ceq = [
        %x - BernsteinAntiDeric(v.*c) + cosntat
        %BernsteinDegreElev(y,2*N-1) - BernsteinAntiDeriv(BernsteinMul(v,s),T,y(1))
        dx - v.*c
        dy - v.*s
        dc + s.*w
        ds - c.*w
        c.^2 + s.^2 - 1
        %v.^2 - dx.^2 - dy.^2
    ];

    c = -v;

end

function [c, ceq] = nonlcon_obstacle(X,init_conds,final_conds,N,T,shapes)
    [~, ceq] = basicnonlcon(X,init_conds,final_conds,N,T);

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
function [X] = InitialGuess(N,T,init_conds,final_conds)

    omg_max = 10*pi/180;
    vmax = norm(final_conds(1:2)-final_conds(1:2))/T*1.5;
    vmax = max([init_conds(5), final_conds(5), vmax]);
    
    x = max(abs([init_conds(1), final_conds(1)]))*rand(N-1,1);
    y = max(abs([init_conds(2), final_conds(2)]))*rand(N-1,1);
    c = rand(N-1,1);
    s = rand(N-1,1);
    v = vmax*ones(N-1,1);
    w = omg_max*2*rand(N-1,1)-omg_max;
    X = [ x;y;c;s;v;w ];

end

% Get Control Points of Trajectories
function [px, py, c, s, dx, dy, dc, ds, v, w] = getControlPoints(X,init_conds,final_conds,N,T)

    persistent old_N old_T d_mat
    if  isempty(old_N) || isempty(old_T) || old_N ~= N || old_T ~= T
        old_N = N;
        old_T = T;
        d_mat = BernsteinDerivElevMat(N,T);
    end

    X = reshape(X(:), [], 6);

    X = [ init_conds; X; final_conds ];
    
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
