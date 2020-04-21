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
fplot(@(t)BernsteinEval(optim_u,T,t),[ 0 T ],'LineWidth',1.3), title('u');
xlabel('Time (s)'), ylabel('Acceleration (m\cdot s^{-2})');
legend('Bernstein Pol.');
subplot(3,1,2), hold on,
fplot(@(t)BernsteinEval(optim_x,T,t),[ 0 T ],'LineWidth',1.3),
scatter(0:T/(order_plus_one-1):T, optim_x), title('x_1');
xlabel('Time (s)'), ylabel('Position (m)');
subplot(3,1,3), hold on,
fplot(@(t)BernsteinEval(optim_v,T,t),[ 0 T ],'LineWidth',1.3), title('x_2');
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


converged_p = fmincon(@(p)final_cost(p),zeros(order_plus_one,2), [],[], Aeq,beq,[],[],rad_dist);
%converged_p = fmincon(@(p)final_cost(p),zeros(order_plus_one,2), [],[], Aeq,beq);




figure
axis equal
hold on

BernsteinPlot(converged_p,T);

for i = 1:size(centre,1)
    scatter(centre(i,1),centre(i,2),'filled');
    circle(centre(i,1),centre(i,2),min_dist(i));
end

legend('curve','control points')