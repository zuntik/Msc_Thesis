% *** Using Bezier curves ***
% Solve minimum time optimal control for Dubin's Car using fmincon
% 
%     Minimize:      J = tf
% 
%     subject to:    dx/dt = v*r1
%                    dy/dt = v*r2
%                    dr1/dt  = -r2*u
%                    dr2/dt  = r1*u
%                    r1^2 + r2^2 = 1
%                    where  R = [r1 -r2 0; r2 r1 0; 0 0 1]; and r1 = cos(psi); r2 = sin(psi);
%                    and dR/dt = R*S([0;0;u])
%

clear; close all;
global N Dm DmSq BN DU
global Nv Nx Nu Nt tf time d0 dF V rMax
global constType

%% USER INPUTS
doUturn = 0;     % if 1, do U-turn boundary conditions
useScaling = 1;  % if 1, scale problem
constType = 2;   % which nl constraint to use for control
                 % 1: difference of analytic squared turn rate
                 % 2: difference of squared analytic turn rate
                 % 3: difference of absolute value of turn rate

%% Note: N >= d0 + dF + 1
d0 = 2;
dF = 2;
N =  50; % number of control points


%% Problem Specifics
Nv = 1; % number of vehicles
Nx = 4; % number of states 
Nu = 1; % number of controls
Nt = N+1; % number of time nodes

%% Constraints
rMax = 0.1;   % rad/s
V = 25;      % m/s constant velocity

%% Boundary Conditions
psi_0 = 0; 
psi_F= 0; 

x0 = [0 0 cos(psi_0) sin(psi_0)]; xF = [0 1000 cos(psi_F) sin(psi_F)];             % position x, y and heading

%% Note: must choose a final time that is feasible!
tf = 2*norm(xF(1:2) - x0(1:2))/V; % 

%% Problem Scaling
if useScaling == 1    
    k = max(abs(xF(1:2)));
    DU = diag([k; k; 1; 1]);  % force xF or yF = 1.0
    TU = 1;                   % don't scale time for min time problem!
    VU = DU/TU;

%     disp(['Scaling Problem with [DU TU VU] = [' num2str([DU TU VU],2) ']...']);
else
    DU = 1;
    TU = 1/rMax;   
    VU = 1;
    k = VU;
    disp('Solving Unscaled Problem...');
end 
    
% scale the problem
tf = tf/TU;
x0 = x0/DU; 
xF = xF/DU;
rMax = rMax * TU;
V = V/k;


time = linspace(0,1,Nt);

%% Bezier Stuff
Dm = Diff_elev(N,1);            % Note: DC = Dm'*C
DmSq = Dm*Dm;                    % Note: DDC = Dm'*Dm'*C
BN = bernsteinMatrix(N,time);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Do the optimization!
% [X, lb, ub] = initGuessBC2(x0, xDot0, xF, xDotF);

% get initial guess and bounds
Cx = zeros(N+1,Nx);
Cx(1,:) = x0;
Cx(end,:) = xF;
for i = 1:N+1
    Cx(i,:) = x0*(N-(i-1))/N + xF*(i-1)/N;
end

Cx = reshape(Cx,[Nx*(N+1),1]);
Cu = zeros(N+1,Nu);

Aeq = zeros(6,(N+1)*(Nx+Nu)+1);

Aeq(1,1) = 1;         %x0
Aeq(2,N+2) = 1;       %y0;
Aeq(3,2*(N+1)+1) = 1; %r10
Aeq(4,3*(N+1)+1) = 1; %r20

Aeq(5,(N+1)) = 1; %xf
Aeq(6,2*(N+1)) = 1; %yf
Aeq(7,3*(N+1)) = 1; %r1f
Aeq(8,4*(N+1)) = 1; %r2f

beq = [x0';xF'];


X = [Cx;Cu;tf]; 

diff = Aeq*X - beq; 

% states
x_lb = -ones(Nx*(N+1),1)*inf;
x_ub = -x_lb;

% control 
u_lb = -ones(Nu*(N+1),1)*rMax;
u_ub = -u_lb; 

% time
t_lb = 1e-4;
t_ub = inf;

lb = [x_lb(:);u_lb;t_lb];
ub = [x_ub(:);u_ub;t_ub];

A = [];
b = [];

options = optimoptions(@fmincon,'Algorithm','sqp',...
                       'MaxFunctionEvaluations',40000,...
                       'ConstraintTolerance',1e-6,...
                       'StepTolerance',1e-6,...
                       'Display','iter');
%                        'OptimalityTolerance',1e-8,...
% %                       'FunctionTolerance',1e-8,..
%                       'ConstraintTolerance',1e-3);

disp('norm of linear contraints before');
disp(norm(Aeq*X-beq));
disp('norm of nonlinear constraints before');
[ c , ceq]= nonlcon(X);
disp(norm(ceq));

tic;
    [xOut,Jout,exitflag,output] = fmincon(@costFunc,X,A,b,Aeq,beq,lb,ub,@nonlcon,options);
toc

disp('norm of linear contraints after');
disp(norm(Aeq*xOut-beq));
disp('norm of nonlinear constraints after');
disp(norm(nonlcon(xOut)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plot results

[XN, YN, r1N, r2N, dXN, dYN, dr1N, dr2N, UN] = getTrajectories(xOut);

% rescale to original values
Temp  = [XN, YN, r1N, r2N]*DU;

XN = Temp(:,1); YN = Temp(:,2);  
r1N = Temp(:,3); r2N = Temp(:,4);
Mag = sqrt(r1N.^2 + r2N.^2);
r1N = r1N./Mag; r2N = r2N./Mag;

tf = Jout*TU;
time = time*TU;
rMax = rMax/TU;
V = V*k;
x0 = x0*DU; 
xF = xF*DU;

Temp  = [dXN, dYN, dr1N, dr2N]*DU;
dXN = Temp(:,1); dYN = Temp(:,2);  
dr1N = Temp(:,3); dr2N = Temp(:,4);
dr1N = dr1N./Mag; dr2N = dr2N./Mag;


BN = bernsteinMatrix(N,time);
XN = BN*XN; 
YN = BN*YN; 
PsiN = BN*atan2(r2N,r1N);

dXN = BN*dXN/tf; 
dYN = BN*dYN/tf;
dPsiN = BN*UN;


% dXN = dXN/tf; dYN = dYN/tf; dPsiN = dPsiN/tf;

time = time*tf;
figure(1);
set(gcf,'Position',[112 123 560 975]);
subplot(6,1,[1:3]);
plot(YN,XN); grid on; axis equal;
title('Path'); ylabel('X (North) [m]');
subplot(6,1,4);
plot(time,PsiN); grid on;
ylabel('Psi [rad]'); %xlabel('Time [s]'); 
subplot(6,1,5);
plot(time,XN); grid on;
ylabel('X (North) [m]'); %xlabel('Time [s]'); 
subplot(6,1,6);
plot(time,YN); grid on;
ylabel('Y (East) [m]'); xlabel('Time [s]'); 

figure(3)
plot(time,dXN); hold
plot(time,dYN);

figure(2);
hLimits = plot(time,ones(size(time))*V,'k--',time,ones(size(time))*rMax,'r:',...
               time,-ones(size(time))*rMax,'r:','LineWidth',3);
grid on; hold on;
hData = plot(time,dXN,time,dYN,time,dPsiN,time,sqrt(dXN.^2 + dYN.^2),...
             'LineWidth',2);
title('Velocity'); xlabel('Time [s]');
hl = legend([hLimits; hData],'$V_{max}$','$r_{max}$','$-r_{max}$',...
            '$\dot{X}$ [m/s]','$\dot{Y}$ [m/s]','r [rad/s]','V  [m/s]');
set(hl, 'Interpreter', 'latex','Location','Best');


save trajectory_final_v1


%% convert control points to state trajectories and their derivatives
function [XN, YN, r1N, r2N, dXN, dYN, dr1N, dr2N, UN] = getTrajectories(X)

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

    global V DU 
    
    tf = X(end);
    
    [~, ~, r1N, r2N, dXN, dYN, dr1N, dr2N, UN] = getTrajectories(X);
  
    % differential constraint
    Mag = sqrt(r1N.^2 + r2N.^2);
    r1N = (r1N + eps)./(Mag+eps); r2N = (r2N + eps)./(Mag + eps);

    ceq = [dXN - tf*V*r1N;dYN - tf*V*r2N; dr1N + tf*r2N.*UN; dr2N - tf*r1N.*UN]; 
%     ceq = [dXN - tf*V*r1N;dYN - tf*V*r2N; dr1N + tf*r2N.*UN; dr2N - tf*r1N.*UN; r1N.^2 + r2N.^2 - 1]; 
    
    
    
%     ceq = [dXN - DU(1,1)\tf*V*cos(PsiN);dYN - DU(2,2)\tf*V*sin(PsiN); dPsiN - DU(3,3)\tf*UN]; 

 
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
function Prod_T = Prod_Matrix(N)
%This function produces a matrix which can be used to compute ||x dot x||^2
% i.e. xaug = x'*x;
% xaug = reshape(xaug',[length(x)^2,1]);
% y = Prod_T*xaug;
% or simply norm_square(x)


T = zeros(2*N+1,(N+1)^2);

for j = 0:2*N
for i = max(0,j-N): min(N,j)
   if N >= i && N >= j-i && 2*N >= j && j-i >= 0
       T(j+1,N*i+j+1) = nchoosek_mod(N,i)*nchoosek_mod(N,j-i)/nchoosek_mod(2*N,j);
   end
end   
end

Prod_T = T;


end


%% Generate n choose k_mod

function out = nchoosek_mod(n,k)

out = 1;
    for i = 1:k
        out = out*(n-(k-i));
        out = out/i;
    end
end
