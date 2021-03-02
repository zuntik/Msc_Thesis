% Written by Thomas Berry
% Inspired by Nguyen T. Hung's code, which in turn is based on the tracking
% controller described in the thesis of Vanni-Vanni 

addpath('..\Bernstein')

% for the medusa
MODELPARAMS.mass = 17.0;
MODELPARAMS.I_z = 1;
MODELPARAMS.X_dot_u = -20;
MODELPARAMS.Y_dot_v = -30;%-1.3175;
MODELPARAMS.N_dot_r = -8.69;% -0.5;
MODELPARAMS.X_u = -0.2;
MODELPARAMS.Y_v = -50;
MODELPARAMS.N_r = -4.14; %-0.1
MODELPARAMS.X_uu = -25;
MODELPARAMS.Y_vv = -0.01;%-101.2776;
MODELPARAMS.N_rr = -6.23; %-21
MODELPARAMS.m_u = MODELPARAMS.mass - MODELPARAMS.X_dot_u;
MODELPARAMS.m_v = MODELPARAMS.mass - MODELPARAMS.Y_dot_v;
MODELPARAMS.m_r = MODELPARAMS.I_z - MODELPARAMS.N_dot_r;
MODELPARAMS.m_uv = MODELPARAMS.m_u - MODELPARAMS.m_v;


% Simulation time
% T = 200;
% cp = [ linspace(0,T,100).', 10.*sin(1.5*2*pi/T*linspace(0,T,100)).'];
T = constants.T;
cp = xOut(:,1:6);

dcp = BernsteinDeriv(cp,T);
ddcp = BernsteinDeriv(dcp,T);
input = [   BernsteinEvaluator(cp(:,1:2), T);
            BernsteinEvaluator(dcp(:,1:2),T);
            BernsteinEvaluator(ddcp(:,1:2),T);
];

% Simulation of the Controlled System

dynamics = @(t,x,dx,u)dynamicspluskinematics_medusa(t,x,dx,u,MODELPARAMS);
controller = @(t,input,x,dx)trajectorytrackingcontroller_medusa(t,input,x,dx,MODELPARAMS);
% dynamics = @dynamicspluskinematics;
% controller = @trajectorytrackingcontroller;

X0 = cp(1,:).';
dX0 = dcp(1,:).';

% X0 = [  input{1}(0); 
%         input{2}(0);
%         atan2(input{4}(0), input{3}(0)); 
%         sqrt(input{4}(0)^2 + input{3}(0)^2);
%         0; 
%         0];
% dX0 = [0;0;0;0;0;0];

% X0 = [ input{1}(0); input{2}(0); atan2(input{4}(0), input{3}(0))];
% dX0 = [0;0;0];
% X0 = [15;-5;pi/2]; % [ px; py; psi];
% dX0 = [0; 0; 0];

[time,Xout] = ode15i(@(t,x,dx) controlled_system(t,x,dx,dynamics,controller,input), [0, T], X0, dX0);


[~,dXout] = gradient(Xout, mean(diff(time)));

Uout=cell(1,length(time));
for i = 1:length(time)
    Uout{i} = controller(time(i),input,Xout(i,:),dXout(i,:));
end
Uout = cell2mat(Uout);


%% Plot
figure, axis equal
hold on;
plot(Xout(:,2),Xout(:,1));
fplot(input{2},input{1},[0 T]);
xlabel('y'),ylabel('x')
if ~isempty(constants.obstacles) && true
    for i = 1:size(constants.obstacles,3)
        plot(polyshape(constants.obstacles(:,[2,1],i)))
    end
end
if ~isempty(constants.obstacles_circles)
    for i = 1:size(constants.obstacles_circles,3)
        centrx = constants.obstacles_circles(i,1);
        centry = constants.obstacles_circles(i,2);
        r = constants.obstacles_circles(i,3);
        plot(polyshape([(centry+r*cos(0:0.01:2*pi));(centrx+r*sin(0:0.01:2*pi))].'));
    end
end
legend('Trajectory Tracking Result','Motion Planning Result');

figure
plot(time,Uout(1,:))
figure
plot(time,Uout(2,:))


% subplot(2,2,1)
% title('\tau_u from Motion Planning')
% subplot(2,2,2)
% title('\tau_r from Motion Planning')
% subplot(2,2,3)
% title('\tau_u from Trajectory Tracking')
% subplot(2,2,4)
% title('\tau_r from Trajectory Tracking')

function F = controlled_system(t,X,dX,dynamics,controller,input)
    % generate the control law
    U = controller(t,input,X,dX);
    % plug it in the dynamics
    F = dX-dynamics(t,X,dX,U);
end

% the system to be controlled
function dX = dynamicspluskinematics(~, X,~, U)
    % a kinematics only model
    % This model is time independent, therefor, the "time" input is ifnored
    ur = U(1);
    r = U(2);
    psi = X(3);
    dX = [  ur*cos(psi); % dx
            ur*sin(psi); % dy
            r];          % dyaw
end

function dX = dynamicspluskinematics_medusa(~, X, ~, U, MODELPARAMS)
    
    % get the values
    yaw = X(3);
    u = X(4);
    v = X(5);
    r = X(6);

    tau_u = U(1);
    tau_r = U(2);

    % load the parameters
    X_u = MODELPARAMS.X_u;
    Y_v = MODELPARAMS.Y_v;
    N_r = MODELPARAMS.N_r;
    X_uu = MODELPARAMS.X_uu;
    Y_vv = MODELPARAMS.Y_vv;
    N_rr = MODELPARAMS.N_rr;

    %%%%%% masses
    m_u = MODELPARAMS.m_u;
    m_v = MODELPARAMS.m_v;
    m_r = MODELPARAMS.m_r;
    m_uv = MODELPARAMS.m_uv;

    %%%%%%%%% drag
    d_u = -X_u - X_uu*abs(u);
    d_v = -Y_v - Y_vv*abs(v);
    d_r = -N_r - N_rr*abs(r);

    dX = zeros(6,1);
    dX(1) = u*cos(yaw) - v*sin(yaw); %dx
    dX(2) = u*sin(yaw) + v*cos(yaw); % dy
    dX(3) = r; % dyaw
    dX(4) = 1/m_u*(tau_u + m_v*v*r - d_u*u); %du
    dX(5) = 1/m_v*(-m_u*u*r - d_v*v); % dv
    dX(6) = 1/m_r*(tau_r + m_uv*u*v - d_r*r); % dr  
    
end

% produces the control law
function U = trajectorytrackingcontroller(t, input, X, ~)
    
    pd = [input{1}(t); input{2}(t)];
    pd_dot = [input{3}(t); input{4}(t)];
    p = X(1:2);
    psi = X(3);
    
    % Setup constraint for the vehicle input (speed and heading rate)
    umin=0;     umax=2; % limit on the vehicle's speed
    rmin=-0.5;  rmax=0.5; % limit on the vehicle's heading rate
    l_bound=[umin;rmin];  u_bound=[umax;rmax];  

    delta=-0.5;
    Delta= [1    0;
            0  -delta];
    epsilon=[delta; 0];        
    kx=.2; ky=0.05;    
    Kk=[kx 0;
        0  ky];
    RB_I=[cos(psi)   -sin(psi);
          sin(psi)    cos(psi)]; 
    e_pos=RB_I'*(p-pd)-epsilon;

    U = Delta\(-Kk*e_pos+RB_I'*pd_dot); 
    U = max(min(U,u_bound),l_bound);

end

function U = trajectorytrackingcontroller_medusa(t,input,X,dX,MODELPARAMS)

    % get the values
    x = X(1);
    y = X(2);
    yaw = X(3);
    u = X(4);
    v = X(5);
    r = X(6);
    p = [x;y];
    psi = yaw;
    pd = [input{1}(t); input{2}(t)];
    pd_dot = [input{3}(t); input{4}(t)];
    pd_dot_dot = [input{5}(t); input{5}(t)];
    
    % outter loop - kinematics controller
    delta=-0.5;
    Delta= [1    0;
            0  -delta];
    epsilon=[delta; 0];        
    kx=.2; ky=0.05;    
    Kk=[kx 0;
        0  ky];
    RB_I=[cos(psi)   -sin(psi);
          sin(psi)    cos(psi)];
    S = [0, -r; r, 0];
    e_pos=RB_I'*(p-pd)-epsilon;
    e_pos_dot = -S*RB_I'*([dX(1);dX(2)]-pd_dot);
    u_d = Delta\(-Kk*tanh(e_pos-[delta;0])-[0;v]+RB_I'*pd_dot); 
    u_d_dot = Delta\( -Kk*e_pos_dot - S*RB_I'*pd_dot + RB_I'*pd_dot_dot );
    
    % inner loop - dynamic controller
    ku = 1;
    kr = 1;
    Kd = [ku, 0; 0, kr];
    
    % bounds
    u_bound = [25.9; .113];
    l_bound = [0; -.113];
    
    % load the parameters
    X_u = MODELPARAMS.X_u;
    N_r = MODELPARAMS.N_r;
    X_uu = MODELPARAMS.X_uu;
    N_rr = MODELPARAMS.N_rr;

    %%%%%% masses
    m_u = MODELPARAMS.m_u;
    m_v = MODELPARAMS.m_v;
    m_r = MODELPARAMS.m_r;
    m_uv = MODELPARAMS.m_uv;

    %%%%%%%%% drag
    d_u = -X_u - X_uu*abs(u);
    d_r = -N_r - N_rr*abs(r);

    %%%%%%%%% Compact Form
    M = [m_u, 0; 0, m_r];
    C = [0, -m_v*v; -m_uv*v, 0];
    D = [d_u, 0; 0, d_r];
    
    U = -Kd*([u;r]-u_d)+M*u_d_dot + C*[u;r] + D*u_d;
    
    U = max(min(U,u_bound),l_bound);
end

