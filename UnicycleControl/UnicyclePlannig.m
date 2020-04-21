addpath('../Bernstein');


% initial conditions
posx0 = 3
posy0 = 5

velx0 = 1
vely0 = 0


% order of polynomial 
n = 10;

% Time 
T = 10;

% the optimisation variable will be a matrix which will be auto flattened
% when calculating the linear constraints 
cp = zeros((n+1) ,4);

cp(1,1:2) = [posx0, posy0];
cp(end,1:2) [0, 0];
cp(2,1:2) = [T*velx0/n+posx0, T*vely0+posy0];


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
BernsteinPlot(p,T,'PlotControlPoints',false);
xlabel('Time (s)');
ylabel('p');
subplot(3,1,2);
BernsteinPlot(v,T,'PlotControlPoints',false);
xlabel('Time (s)');
ylabel('v');
subplot(3,1,3);
BernsteinPlot(u,T,'PlotControlPoints',false);
xlabel('Time (s)');
ylabel('u');
sgtitle({
    ['unicycle']
    ['Computation time = ' num2str(elapsed_time) ' seconds']});
savefig('FullyDefined.fig');
