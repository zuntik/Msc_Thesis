clear all
close all
addpath('../BeBOT_lib');
addpath('../Bernstein');

%% Problem

pos0 = 3;
vel0 = 1;
init_cond = [ pos0; vel0];

T = 10;

%oder of polynomial
n = 10;

%% Optimize tracked trajectories
% i.e. integrate trajectories obtained via trajectory tracking


% cost function for x which needs trajectory tracking
% the only constraints are the inital state constraints, which can be
% the dynamic contraints are imposed right when the cost is calculated.
% ode45 is used to simulate the remaining trajectories implemented as
% linear constraints

% gains for the feedback loop
Kp = 1;
Kv = 10;


cp = zeros(n+1,1);
cp(1) = pos0;
cp(2) = T*vel0/n + pos0;
Aeq_tt = zeros(3,length(cp));
Aeq_tt(1,1) = 1;
Aeq_tt(2,1) = -n/T;
Aeq_tt(2,2) = n/T;
Aeq_tt(3,end)=1;
beq_tt = [ pos0 ; vel0; 0];

tic
cp = fmincon(@(x)costfun_tt(x,T,Kp,Kv), cp, [], [], Aeq_tt, beq_tt);
elapsed_time = toc;


dcp  = BernsteinDeriv(cp, T);
ddcp = BernsteinDeriv(dcp,T);
pd = @(t) BernsteinEval(cp,  T,t);
vd = @(t) BernsteinEval(dcp, T,t);
ad = @(t) BernsteinEval(ddcp,T,t);

[t,y] = ode45(@(t,y) TrajectoryTrackingODE(t,y,pd,vd,ad,Kp,Kv), [0 T], [cp(1);(cp(2)-cp(1))*n/T;0]);
p = y(:,1);
v = y(:,2);
pd = arrayfun(pd,t);
vd = arrayfun(vd,t);
ad = arrayfun(ad,t);
vref = vd - Kp.*(p-pd);
u = v.*abs(v) - Kv.*(v - vref) + ad - Kp.*(v-vd);

figure
subplot(3,1,1);
plot(t,p);
xlabel('Time (s)');
ylabel('p');
subplot(3,1,2);
plot(t,v);
xlabel('Time (s)');
ylabel('v');
subplot(3,1,3);
plot(t,u);
xlabel('Time (s)');
ylabel('u');


sgtitle({
    ['Trajectory Tracking method']
    ['Computation time = ' num2str(elapsed_time) ' seconds']});


%% Optimize planned trajectories
% fully defined variables
% i.e. assume variables are given by polynomials and perform calculations
% on them


% the initial and final 

cp = zeros((n+1)*3 ,1);
cp(1) = pos0;
cp(2) = T*vel0/n + pos0;
cp(n+2:2*n+2) = BernsteinDerivMat(n+1,T)*BernsteinDegrElevMat(n,n+1)*cp(1:n+1);

[Aeq_fd, beq_fd] = LinearConstr_fd(n,T,pos0,vel0);

cp(n+2:2*n+2) = BernsteinDerivMat(n+1,T)*BernsteinDegrElevMat(n,n+1)*cp(1:n+1);

tic 
options = optimoptions('fmincon','Algorithm','sqp');
cp = fmincon(@(x)costfun_fd(x,n,T), cp, [], [], Aeq_fd, beq_fd, [],[], @(x)nonlcon_fd(x,n,T,1,1),options);
elapsed_time = toc;

p = cp(1:n+1);
v = cp(n+2: 2*n+2);
u = cp(2*n+3: 3*n+3);
figure, hold on,
subplot(3,1,1);
BernsteinPlot(p,T,'PlotControlPoints',true);
xlabel('Time (s)');
ylabel('p');
subplot(3,1,2);
BernsteinPlot(v,T,'PlotControlPoints',true);
xlabel('Time (s)');
ylabel('v');
subplot(3,1,3);
BernsteinPlot(u,T,'PlotControlPoints',true);
xlabel('Time (s)');
ylabel('u');
sgtitle({
    ['Point Evaluation Method']
    ['Computation time = ' num2str(elapsed_time) ' seconds']});
savefig('FullyDefined.fig');
