% *** Using Bezier curves ***
% Solve minimum time optimal control for Dubin's Car using fmincon
% 
%     Minimize:      J = rotation speed + linear acceletaion 
% 
%     subject to:    dx/dt = v*r1
%                    dy/dt = v*r2
%                    dv/dt = u2
%                    dr1/dt  = -r2*u1
%                    dr2/dt  = r1*u1
%                    r1^2 + r2^2 = 1
%                    where  R = [r1 -r2 0; r2 r1 0; 0 0 1]; and r1 = cos(psi); r2 = sin(psi);
%                    and dR/dt = R*S([0;0;u1])

clear; close all;

addpath('Bernstein')

global N Nt 
global Nv Nx Nu T time rMax vMax

%% USER INPUTS
doUturn    = 0; % if 1, do U-turn boundary conditions
useScaling = 0; % if 1, scale problem

N =  50; % order


%% Problem Specifics
Nv = 1;   % number of vehicles
Nx = 5;   % number of states 
Nu = 2;   % number of controls
Nt = N+1; % number of time nodes (= n� of control points)

%% Constraints
rMax = 0.1;   % rad/s
vMax = 1;   % m/s

%% Boundary Conditions
psi_0 = 0; 
psi_F = 0; 

% position x, y and heading
x0 = [0 0 cos(psi_0) sin(psi_0)];  
xF = [0 5 cos(psi_F) sin(psi_F)];
v0 = 0.1;
vF = 1;

%% Note: must choose a final time that is feasible!
T = 2*norm(xF(1:2) - x0(1:2))/vMax; % 


time = linspace(0,1,Nt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Do the optimization!
% [X, lb, ub] = initGuessBC2(x0, xDot0, xF, xDotF);

% get initial guess and bounds
Cx = zeros(N+1,Nx);
Cx(1,:) = [x0, v0];
Cx(end,:) = [xF, vF];
for i = 1:N+1
    Cx(i,:) = [x0, v0]*(N-(i-1))/N + [xF, vF]*(i-1)/N;
end

Cx = reshape(Cx,[Nx*(N+1),1]);
Cu = zeros((N+1)*Nu,1);

Aeq_boundary = zeros(8,(N+1)*(Nx+Nu));

Aeq_boundary(1,1)         = 1; % x0
Aeq_boundary(2,1*(N+1)+1) = 1; % y0
Aeq_boundary(3,2*(N+1)+1) = 1; % r10
Aeq_boundary(4,3*(N+1)+1) = 1; % r20

Aeq_boundary(5,1*(N+1))   = 1; % xf
Aeq_boundary(6,2*(N+1))   = 1; % yf
Aeq_boundary(7,3*(N+1))   = 1; % r1f
Aeq_boundary(8,4*(N+1)+1) = 1; % r2f

Aeq_dyn = [ zeros(N+1,(N+1)*4), BernsteinDerivElevMat(N,T), zeros(N+1), -eye(N+1) ]; % dv = u2
beq_dyn = zeros(N+1,1);

Aeq = [Aeq_boundary; Aeq_dyn];
beq = [x0';xF';beq_dyn];

size(Aeq_boundary)
size(Aeq_dyn)
size(Aeq)
size(beq)

X = [Cx;Cu]; 

% states
x_ub = ones(Nx*(N+1),1)*inf;
x_lb = -x_ub;

% control 
u_ub = ones(Nu*(N+1),1)*rMax;
u_lb = -u_ub; 

lb = [x_lb(:); u_lb(:)];
ub = [x_ub(:); u_ub(:)];

A = [];
b = [];

options = optimoptions(@fmincon,'Algorithm','sqp',...
                        'MaxFunctionEvaluations',40000,...
                        'ConstraintTolerance',1e-6,...
                        'StepTolerance',1e-6,...
                        'Display','iter');
%                        'OptimalityTolerance',1e-8,...
%                        'FunctionTolerance',1e-8,..
%                        'ConstraintTolerance',1e-3);

tic;
    [xOut,Jout,exitflag,output] = fmincon(@costFunc,X,A,b,Aeq,beq,lb,ub,@nonlcon,options);
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plot results

[XN, YN, r1N, r2N, VN, dXN, dYN, dr1N, dr2N, dVN, UN1, UN2] = getTrajectories(xOut);

%BN = bernsteinMatrix(N,time);
%XN = BN*XN; 
%YN = BN*YN; 
%PsiN = BN*atan2(r2N,r1N);

%dXN   = BN*dXN/tf; 
%dYN   = BN*dYN/tf;
%dPsiN = BN*UN;


% dXN = dXN/tf; dYN = dYN/tf; dPsiN = dPsiN/tf;

%time = time*tf;
figure(1);
set(gcf,'Position',[112 123 560 975]);
subplot(6,1,1:3);
%plot(YN,XN); grid on; axis equal;
BernsteinPlot([XN, YN],T);
title('Path'); ylabel('X (North) [m]');
%subplot(6,1,4);
%plot(time,PsiN); grid on;
%ylabel('Psi [rad]'); %xlabel('Time [s]'); 
subplot(6,1,5);
plot(time,XN); grid on;
ylabel('X (North) [m]'); %xlabel('Time [s]'); 
subplot(6,1,6);
plot(time,YN); grid on;
ylabel('Y (East) [m]'); xlabel('Time [s]'); 

%figure(3)
%plot(time,dXN); hold
%plot(time,dYN);
%
%figure(2);
%hLimits = plot(time,ones(size(time))*V,'k--',time,ones(size(time))*rMax,'r:',...
%               time,-ones(size(time))*rMax,'r:','LineWidth',3);
%grid on; hold on;
%hData = plot(time,dXN,time,dYN,time,dPsiN,time,sqrt(dXN.^2 + dYN.^2),...
%             'LineWidth',2);
%title('Velocity'); xlabel('Time [s]');
%hl = legend([hLimits; hData],'$V_{max}$','$r_{max}$','$-r_{max}$',...
%            '$\dot{X}$ [m/s]','$\dot{Y}$ [m/s]','r [rad/s]','V  [m/s]');
%set(hl, 'Interpreter', 'latex','Location','Best');


save trajectory_final_v1


%% convert control points to state trajectories and their derivatives
function [XN, YN, r1N, r2N, VN, dXN, dYN, dr1N, dr2N, dVN, UN1, UN2] = getTrajectories(X)

    global N Nx Nu Dm
          
    % get control points
    Cx = X(1:Nx*(N+1),1);
    Cx = reshape(Cx,[N+1,Nx]);

    Cu = X(Nx*(N+1)+1:end-1);
    UN = reshape(Cu,[N+1,Nu]);   
    
    % get states
    XN = Cx(:,1);
    YN = Cx(:,2);
    r1N = Cx(:,3);
    r2N = Cx(:,4);
    VN = Cx(:,5);

    
    % compute  derivatives 
    DCx = Dm'*Cx;
    dXN = DCx(:,1);
    dYN = DCx(:,2);
    dr1N = DCx(:,3);
    dr2N = DCx(:,4);

    
end


%% Nonlinear constraints
% NOTE: should I penalize the entire trajectory or just the max values? 
function [c, ceq] = nonlcon(X)

    [~, ~, r1N, r2N, VN, dXN, dYN, dr1N, dr2N, dVN, UN1, UN2] = getTrajectories(X);
  
    % differential constraint
    Mag = sqrt(r1N.^2 + r2N.^2);
    
    ceq = [ 
        dXN - V.*r1N
        dYN - V.*r2N
        dr1N + r2N.*UN1
        dr2N - r1N.*UN1
        Mag - 1
    ]; 
 
    c = []; 
 
end


%% Cost function
function J = costFunc(X)    
    J = X(end);
end    



%% Create initial guess and lower/upper bound vectors
 
% Bezier derivative
function Dm = Diff(N,tf )
% derivative of a Bezier curve
% INPUT
% N: number of nodes
% OUTPUT
% Dm{N}: differentiation matrix for bez curves of order N (N+1 ctrl points)

% Notes:
% If Cp are the control points of bez, then the control points of bezdot are Cpdot = Cp*Dm
% To compute bezdot with N+1 ctrl points, degree elevation must be performed 

    Dm = -[N/tf*eye(N); zeros(1,N)]+[zeros(1,N);N/tf*eye(N)];

end

function Telev = deg_elev(N)
% INPUT
% N order of Bezier curve
% OUTPUT
% Telev{N}: Transformation matrix from Nth order (N+1 control points) to
% (N+1)th order (N+2 control points)
% If Cp is of order N-1, then Cp*Telev{N-1} is of order N
% see Equation (12+1) in https://pdfs.semanticscholar.org/f4a2/b9def119bd9524d3e8ebd3c52865806dab5a.pdf
% Paper: A simple matrix form for degree reduction of Be�zier curves using Chebyshev�Bernstein basis transformations

    if N < 5 
      es='ERROR: The approximation order should be at least 5';
      disp(es); Dm = [];
    return
    end

    for i = 1:1:N
        Telev{i} = zeros(i+2,i+1);
        for j = 1:1:i+1
            Telev{i}(j,j) = i+1-(j-1);
            Telev{i}(j+1,j) = 1+(j-1);
        end
        Telev{i} = 1/(i+1)*Telev{i}';
    end

end




function Dm = Diff_elev(N,tf)
% derivative of a Bezier curve
% INPUT
% N: number of nodes, tf: final time
% OUTPUT
% Dm{N}: differentiation matrix for bez curves of order N (N+1 ctrl points)
% The differentiation matrix is (N+1)by(N+1), ie. differently from Diff,
% this matrix gives a derivative of the same order of the curve to be
% differentiated

% Notes:
% If Cp are the control points of bez, then the control points of bezdot are Cpdot = Cp*Dm

    Dm = Diff(N,tf);
    Telev = deg_elev(N);
    Dm = Dm*Telev{N-1};

end


%% Generate Product Matrix
function Mat = Prod_Matrix(N)
%This function produces a matrix which can be used to compute ||x dot x||^2
% i.e. xaug = x'*x;
% xaug = reshape(xaug',[length(x)^2,1]);
% y = Mat*xaug;
% or simply norm_square(x)

    Mat = zeros(2*N+1,(N+1)^2);

    for j = 0:2*N
        for i = max(0,j-N): min(N,j)
           if N >= i && N >= j-i && 2*N >= j && j-i >= 0
               Mat(j+1,N*i+j+1) = nchoosek_mod(N,i)*nchoosek_mod(N,j-i)/nchoosek_mod(2*N,j);
           end
        end
    end

end


%% Generate n choose k_mod

function out = nchoosek_mod(n,k)
    out = 1;
    for i = 1:k
        out = out*(n-(k-i));
        out = out/i;
    end
end
