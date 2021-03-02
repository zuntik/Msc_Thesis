% Written by Thomas Berry
% Inspired by Nguyen T. Hung's code, which in turn is based on the tracking
% controller described in the thesis of Vanni-Vanni 

% Simulation time
T = 200;

% Function of time for the input (the desired position and it's derivative)
input =  {  @(t) t;  % x position
            @(t) 10.*sin(0.05*t); % y position
            @(t) 1; % dx
            @(t) 0.5.*cos(0.05*t);};  % dy

% Simulation of the Controlled System

X0 = [15;-5;pi/2]; % [ x; y; psi];

[time,Xout] = ode45(@(t,X) controlled_system(t,X,@dynamics,@controller,input), [0, T], X0);


%% Plot
figure
hold on;
plot(Xout(:,1),Xout(:,2));
fplot(input{1},input{2},[0 T]);
legend('vehicle trajectory', 'desired trajectory'); 

function dX = controlled_system(t,X,dynamics,controller,input)
    % generate the control law
    U = controller(t,input,X);
    % plug it in the dynamics
    dX = dynamics(t,X,U);
end

% the system to be controlled
function dX = dynamics(~, X, U)
    % a kinematics only model
    % This model is time independent, therefor, the "time" input is ifnored
    ur = U(1);
    r = U(2);
    psi = X(3);
    dX = [  ur*cos(psi); % dx
            ur*sin(psi); % dy
            r];          % dyaw
end

% produces the control law
function U = controller(t, input, X)
    
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
