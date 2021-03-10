clear all
close all

    
%% Load Parameters    
CONSTANTS.N = 20; % Order of approximation
CONSTANTS.pinit = [0;0];
CONSTANTS.pfin = [10;10];
CONSTANTS.pobs = [5;5];
CONSTANTS.psiinit = pi/4;
CONSTANTS.psifin = pi/4;
CONSTANTS.vinit = 0.1;
CONSTANTS.vfin = 0.1;
CONSTANTS.omegainit = 0;
CONSTANTS.omegafin = 0;
CONSTANTS.vmax = 5;
CONSTANTS.sep = 0.5;
CONSTANTS.omegamax = 1;
CONSTANTS.T =100;


%% Initial Guess
x = linspace(CONSTANTS.pinit(1),CONSTANTS.pfin(1),CONSTANTS.N+1)'+rand;
y = linspace(CONSTANTS.pinit(1),CONSTANTS.pfin(1),CONSTANTS.N+1)';
psi = linspace(CONSTANTS.psiinit(1),CONSTANTS.psifin(1),CONSTANTS.N+1)';
omega = linspace(CONSTANTS.omegainit(1),CONSTANTS.omegafin(1),CONSTANTS.N+1)';
v = linspace(CONSTANTS.vinit(1),CONSTANTS.vfin(1),CONSTANTS.N+1)';

x0 = [x;y;psi;omega;v];


%% Linear Constraints
N = CONSTANTS.N;
A=[]; b=[]; Aeq=[]; beq=[]; 
lb=[
        -Inf*ones(N+1,1);
        -Inf*ones(N+1,1);
        -Inf*ones(N+1,1);
        -CONSTANTS.omegamax*ones(N+1,1);
        0*ones(N+1,1);
    ]; 
 
ub=[
        Inf*ones(N+1,1);
        Inf*ones(N+1,1);
        Inf*ones(N+1,1);
        CONSTANTS.omegamax*ones(N+1,1);
        CONSTANTS.vmax*ones(N+1,1);
    ]; 


%% Solve problem 
% Method 1 for Bernstein, Method 0 for LGL PS
CONSTANTS.method = 1;

options = optimoptions(@fmincon, ...
    'MaxIter', 10000, ...     Maximum number of iterations
    'MaxFunEvals', 100000, ...   Maximum number of function evaluations
    'Display', 'iter-detailed', ... Display the progress
    'DiffMinChange', 1e-5, ...  Similar to the relaxation bound of the problem
    'Algorithm', 'SQP' ...  We are using the SQP algorithm to solve this OCP
    );

tic
[xbern,fbern] = fmincon(@(xbern)costfun(xbern,CONSTANTS),x0,A,b,Aeq,beq,lb,ub,@(xbern)nonlcon(xbern,CONSTANTS),options);
toc

fbern


















%% Plot Results
N = CONSTANTS.N;
T = CONSTANTS.T;

x = xbern(1:N+1);
y = xbern(N+2:2*N+2);
psi = xbern(2*N+3:3*N+3);
omega = xbern(3*N+4:4*N+4);
v = xbern(4*N+5:5*N+5);
[tnodes,w,Diff] = BeBOT(N,T);
time = linspace(0,T,1000);



figure(1)
subplot(2,2,1)
plot(x,y,'o','Color','b','Linewidth',1); hold on
plot(BernsteinPoly(x',time),BernsteinPoly(y',time),'-','Color','b','Linewidth',1);
pos = [
        CONSTANTS.pobs(1)-CONSTANTS.sep ...
        CONSTANTS.pobs(2)-CONSTANTS.sep ...
        2*CONSTANTS.sep 2*CONSTANTS.sep
       ]; 
rectangle('Position',pos,'Curvature',[1 1])
axis equal
title('Bernstein method - Trajectory')
grid on

subplot(2,2,2)
plot(tnodes,psi,'o','Color','b','Linewidth',1); hold on
plot(time,BernsteinPoly(psi',time),'-','Color','b','Linewidth',1);
title('Bernstein method - Heading')
grid on

subplot(2,2,3)
plot(tnodes,omega,'o','Color','b','Linewidth',1); hold on
plot(time,BernsteinPoly(omega',time),'-','Color','b','Linewidth',1);
title('Bernstein method - Ang Rate')
grid on

subplot(2,2,4)
plot(tnodes,v,'o','Color','b','Linewidth',1); hold on
plot(time,BernsteinPoly(v',time),'-','Color','b','Linewidth',1);
title('Bernstein method - Speed')
grid on































%% Cost Function

function J = costfun(xbern,CONSTANTS)
%COSTFUN Summary of this function goes here
N = CONSTANTS.N; 
T = CONSTANTS.T;

x = xbern(1:N+1);
y = xbern(N+2:2*N+2);
psi = xbern(2*N+3:3*N+3);
omega = xbern(3*N+4:4*N+4);
v = xbern(4*N+5:5*N+5);

if CONSTANTS.method == 1
    [~,w,Diff] = BeBOT(N,T);
end 

if CONSTANTS.method == 0
    [~,w,Diff] = LGL_PS(N,T);
end 

domega = Diff'*omega;
dv = Diff'*v;

J = w'*(2*domega.^2+dv.^2);
end


%% Nonlinear Constraints
function [c,ceq] = nonlcon(xbern,CONSTANTS)

N = CONSTANTS.N; 
T = CONSTANTS.T;

x = xbern(1:N+1);
y = xbern(N+2:2*N+2);
psi = xbern(2*N+3:3*N+3);
omega = xbern(3*N+4:4*N+4);
v = xbern(4*N+5:5*N+5);

if CONSTANTS.method == 1
    [~,~,Diff] = BeBOT(N,T);
end 

if CONSTANTS.method == 0
    [~,~,Diff] = LGL_PS(N,T);
end 

dist2obs2 = (BernsteinProduct(x'-CONSTANTS.pobs(1),x'-CONSTANTS.pobs(1))+BernsteinProduct(y'-CONSTANTS.pobs(2),y'-CONSTANTS.pobs(2)));

c=[
        -dist2obs2'+CONSTANTS.sep^2;
  ];
ceq=[
        x(1)-CONSTANTS.pinit(1);y(1)-CONSTANTS.pinit(2);
        x(end)-CONSTANTS.pfin(1);y(end)-CONSTANTS.pfin(2);
        psi(1)-CONSTANTS.psiinit;
        psi(end)-CONSTANTS.psifin;
        v(1)-CONSTANTS.vinit;
        v(end)-CONSTANTS.vfin;
        omega(1)-CONSTANTS.omegainit;
        omega(end)-CONSTANTS.omegafin;
        Diff'*x - v.*cos(psi);
        Diff'*y - v.*sin(psi);
        Diff'*psi - omega;
    ];

end


