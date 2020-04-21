%%% IIEEC - Thomas Berry


%% 1D problem (a double integrator):


initial_position = 3;
initial_velocity = 1;
x0_1d = [ initial_position ; initial_velocity ];

rho = 0.1;


T = 6;

A = [ 
    0 1
    0 0 
    ];
B = [ 0; 1];
C = [ 1 0 ];
D = 0;
G = ss(A,B,C,D);


%% analytical soltuion to 1D problem



Q = eye(2);
K = lqr(G,Q,rho);

% sys is the closed loop system 
closed_G  = ss(A-B*K, zeros(size(A)),C,D);

[~,t_lqr,x_lqr] = initial(closed_G,x0_1d, T);
u_lqr = -K*x_lqr';


%sgtitle('Response of a LQR controlled double Integrator');
subplot(3,1,1), plot(t_lqr,u_lqr,'Linewidth',1.3), title('u');
xlabel('Time (s)'), ylabel('Acceleration (m\cdot s^{-2})');
subplot(3,1,2), plot(t_lqr,x_lqr(:,1),'Linewidth',1.3), title('x_1');
xlabel('Time (s)'), ylabel('Position (m)');
subplot(3,1,3), plot(t_lqr,x_lqr(:,2),'Linewidth',1.3), title('x_2');
xlabel('Time (s)'), ylabel('Velocity (m\cdot s^{-1})');


%% zoh
T = 6;
% choose order for u
order = 30;
coefs = rand(1,order);
% define cost
Th = T/order;
cost=@(coefs) cost_zoh(coefs, Th, initial_velocity, initial_position, rho);
    
% search 
options=optimset('MaxFunEvals',1e100*length(coefs), 'MaxIter',1e100*length(coefs));

tic
for i = 1:5
    coefs = fminsearch(cost,coefs,options);
end
disp(toc)
    
% plot solution
Ts = 0.01;
t = 0:Ts:T-Ts;
u = repelem(coefs,Th/Ts);
[~,t,x] = lsim(G,u,t,x0_1d,'zoh');

%%

subplot(3,1,1), hold on,
plot(t_lqr,u_lqr,'Linewidth',1.3), plot(t,u,'Linewidth',1.3), title('u');
legend('analytic','shooting');
xlabel('Time (s)'), ylabel('Acceleration (m\cdot s^{-2})');
subplot(3,1,2), hold on,
plot(t_lqr,x_lqr(:,1),'LineWidth',1.3), plot(t,x(:,1),'Linewidth',1.3), title('x_1');
xlabel('Time (s)'), ylabel('Position (m)');
subplot(3,1,3), hold on,
plot(t_lqr,x_lqr(:,2),'LineWidth',1.3), plot(t,x(:,2),'Linewidth',1.3), title('x_2');
xlabel('Time (s)'), ylabel('Velocity (m\cdot s^{-1})');


%% poly 

initial_position = 3;
initial_velocity = 1;
x0_1d = [ initial_position ; initial_velocity ];

rho = 0.1;


% choose order for pol of x
T = 6;
order_plus_one =8;

% search
Aeq = zeros(2,order_plus_one);
Aeq(1,end) = 1;
Aeq(2,end-1) = 1;
beq = x0_1d;
Aeq(3,:) = T.^[order_plus_one-1:-1:0];
beq(3) = 0;

optim_x = [x0_1d', rand(1,order_plus_one-2)];

tic
for i = 1:5
    optim_x = fmincon(@(p)(cost_pol2(p,T,rho)), optim_x, [],[], Aeq, beq);
end
disp(toc);
optim_v = polyder(optim_x);
optim_u = polyder(optim_v);


subplot(3,1,1), hold on,
plot(t_lqr,u_lqr,'Linewidth',1.3),  fplot(@(t)polyval(optim_u,t),[ 0 T ],'LineWidth',1.3), title('u');
xlabel('Time (s)'), ylabel('Acceleration (m\cdot s^{-2})');
legend('analytic','Polynomial')
subplot(3,1,2), hold on,
plot(t_lqr,x_lqr(:,1),'LineWidth',1.3), fplot(@(t)polyval(optim_x,t),[ 0 T ],'LineWidth',1.3), title('x_1');
xlabel('Time (s)'), ylabel('Position (m)');
subplot(3,1,3), hold on,
plot(t_lqr,x_lqr(:,2),'LineWidth',1.3), fplot(@(t)polyval(optim_v,t),[ 0 T ],'LineWidth',1.3), title('x_2');
xlabel('Time (s)'), ylabel('Velocity (m\cdot s^{-1})');


%% Linear quadratic Programming

Umax = 1;

% simulation time
T = 6;

% sample time
Ts = 0.1;

%number of samples
N = T/Ts;

t = 0:Ts:T;

rho = .1;
Q = eye(2);

% Descritized G
dG = c2d(G,Ts);
Phi = dG.A;
Gamma = dG.B;

dim = size(Phi,1);
dim_u = size(Gamma,2);
dim_xbar = N*(dim+dim_u) + dim;

AeqDynamics = [];

H = [];
QandR = blkdiag(Q,rho);

Phi_Gamma_plus1 = [ Phi, Gamma, -eye(dim) ];

for i = 1:N
    line_block = [zeros(dim, (dim+dim_u)*(i-1)) , Phi_Gamma_plus1, zeros(dim, (dim+dim_u)*(N-i)) ];
    AeqDynamics = [ AeqDynamics ;line_block];
    H = blkdiag(H,QandR);
end



H = blkdiag(H,Q);
H= 1/2 * H;

AeqInitialState = [eye(dim), zeros(dim, N*(dim+dim_u))];
AeqFinalState   = [zeros(dim, N*(dim+dim_u)),eye(dim)];
Aeq = [AeqDynamics; AeqInitialState; AeqFinalState];

ub = Umax  * [repmat([Inf*ones(1,dim),1], 1,N), Inf*ones(1,dim)];
lb = -Umax * [repmat([Inf*ones(1,dim),1], 1,N), Inf*ones(1,dim)];

beq = [ zeros(N*dim,1) ; x0_1d ; zeros(dim,1) ];

f = zeros(1,dim_xbar);

tic
x_bar_opt = quadprog(H,f,[],[],Aeq,beq,lb,ub);
toc

clear x
x(:,1) = x_bar_opt(1:3:end);
x(:,2) = x_bar_opt(2:3:end);
u =   x_bar_opt(3:3:end);

figure
if Umax == Inf
    subplot(3,1,1), hold on,
    plot(t_lqr,u_lqr,'Linewidth',1.3), plot(t(1:end-1),u,'Linewidth',1.3), title('u');
    xlabel('Time (s)'), ylabel('Acceleration (m\cdot s^{-2})');
    legend('analytic','Quad. Prog');
    subplot(3,1,2), hold on,
    plot(t_lqr,x_lqr(:,1),'LineWidth',1.3), plot(t,x(:,1),'Linewidth',1.3), title('x_1');
    xlabel('Time (s)'), ylabel('Position (m)');
    subplot(3,1,3), hold on,
    plot(t_lqr,x_lqr(:,2),'LineWidth',1.3), plot(t,x(:,2),'Linewidth',1.3), title('x_2');
    xlabel('Time (s)'), ylabel('Velocity (m\cdot s^{-1})');
else
    subplot(3,1,1), plot(t(1:end-1),u,'Linewidth',1.3), title('u');
    xlabel('Time (s)'), ylabel('Acceleration (m\cdot s^{-2})');
    subplot(3,1,2), plot(t,x(:,1),'Linewidth',1.3), title('x_1');
    xlabel('Time (s)'), ylabel('Position (m)');
    subplot(3,1,3), plot(t,x(:,2),'Linewidth',1.3), title('x_2');
    xlabel('Time (s)'), ylabel('Velocity (m\cdot s^{-1})');
end

%% Bernstein Polyinomial


addpath('Bernstein');

T = 6;
rho = 0.1;

initial_position = 3;
initial_velocity = 1;
x0_1d = [ initial_position ; initial_velocity ];

% choose order for pol of x
order_plus_one = 15;

Aeq= zeros(2,order_plus_one);
Aeq(1,1) = 1;
Aeq(2,1) = -(order_plus_one-1)/T;
Aeq(2,2) = (order_plus_one-1)/T;
beq = x0_1d;

Aeq(3,end)=1;
beq(3) = 0;

b_sq_int = @(p,T) BernsteinIntegr(BernsteinPow(p,2),T);

cost = @(p) b_sq_int(p,T) + b_sq_int(BernsteinDeriv(p,T),T) + rho* b_sq_int(BernsteinDeriv(BernsteinDeriv(p,T),T),T);

optim_x = zeros(order_plus_one,1);

tic
for i = 1:1
    optim_x = fmincon(cost,optim_x, [],[], Aeq,beq);
end
disp(toc)

optim_v = BernsteinDeriv(optim_x,T);
optim_u = BernsteinDeriv(optim_v,T);

figure
subplot(3,1,1), hold on,
plot(t_lqr,u_lqr,'Linewidth',1.3), fplot(@(t)BernsteinEval(optim_u,T,t),[ 0 T ],'LineWidth',1.3), title('u');
xlabel('Time (s)'), ylabel('Acceleration (m\cdot s^{-2})');
legend('analytic','Bernstein Pol.');
subplot(3,1,2), hold on,
plot(t_lqr,x_lqr(:,1),'LineWidth',1.3), fplot(@(t)BernsteinEval(optim_x,T,t),[ 0 T ],'LineWidth',1.3),
scatter(0:T/(order_plus_one-1):T, optim_x), title('x_1');
xlabel('Time (s)'), ylabel('Position (m)');
subplot(3,1,3), hold on,
plot(t_lqr,x_lqr(:,2),'LineWidth',1.3), fplot(@(t)BernsteinEval(optim_v,T,t),[ 0 T ],'LineWidth',1.3), title('x_2');
xlabel('Time (s)'), ylabel('Velocity (m\cdot s^{-1})');


%% 2D problem (2 double integrators):


initial_positionx = 3;
initial_velocityx = 1;
initial_positiony = 5;
initial_velocityy = 0;

x0_2d = [ 
    initial_positionx initial_positiony
    initial_velocityx initial_velocityy 
];


%% Bernstein Polyinomial 2D
close all
addpath('Bernstein');

order_plus_one = 9;
T=10;

Aeq_1d = zeros(3,order_plus_one);
Aeq_1d(1,1) = 1;
Aeq_1d(2,1) = -(order_plus_one-1)/T;
Aeq_1d(2,2) = (order_plus_one-1)/T;
% Aeq_1d(2,1) = -1;
% Aeq_1d(2,2) = 1;
Aeq_1d(3,end)=1;
Aeq = blkdiag(Aeq_1d,Aeq_1d);
beq = [ initial_positionx; initial_velocityx; 0;...
        initial_positiony; initial_velocityy; 0 ];

rho = 100;

centre = [ 1 1; 3 2.5];
min_dist = [1 ; 0.5];

rad_dist = @(p) rad_dist_const(p,T,centre,min_dist);

b_sq_int = @(p,T) BernsteinIntegr(BernsteinMul(p,p),T);
twoD_cost = @(p,T) b_sq_int(p,T) + b_sq_int(BernsteinDeriv(p,T),T) + rho* b_sq_int(BernsteinDeriv(BernsteinDeriv(p,T),T),T);
final_cost = @(p) norm(twoD_cost(p,T));

tic
converged_p = fmincon(@(p)final_cost(p),zeros(order_plus_one,2), [],[], Aeq,beq,[],[],rad_dist);
%converged_p = fmincon(@(p)final_cost(p),zeros(order_plus_one,2), [],[], Aeq,beq);
toc
%%
figure
%axis([-1 6 ],'DataAspectRatio',[1 1 1])
axis equal
hold on

BernsteinPlot(converged_p,T);

for i = 1:size(centre,1)
    scatter(centre(i,1),centre(i,2),'filled');
    circle(centre(i,1),centre(i,2),min_dist(i));
end

legend('curve','control points')
