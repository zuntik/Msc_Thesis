
%% configs
Kp = 1;
Kv = 10;
Kg = 1;

%% example traj

pol_velocity = [ -1 4 0];
pol_position = [ -1/3 2 0 0 ];
pol_acceleration =  [ -2 4 ];
T = 4;

y0 = [ -1 ; 0];


%% code version

pd = @(t)polyval(pol_position,t);
vd = @(t)polyval(pol_velocity,t);
ad = @(t)polyval(pol_acceleration,t);

[t,y] = ode45(@(t,y) odefcn(t,y,pd,vd,ad,Kp,Kv), [0 4], y0);

figure, hold on, title('position');
plot(t,y(:,1));
fplot(pd, [0 4]);
figure, hold on, title('speed');
plot(t,y(:,2));
fplot(vd,[0 4]);


%% seperatly
%position

[t,y] = ode45(@(t,y) vd(t)- Kp*(y-pd(t)) , [0 4], -1);
figure, hold on
plot(t,y);
fplot(pd,[0 4]);

%speed

[t,y] = ode45(@(t,y) ad(t) - Kv * (y-vd(t)), [0 4],vd(0));
figure, hold on
plot(t,y);
fplot(vd,[0 4]);

%% function for joint

function dydt = odefcn(t,y,pd,vd,ad,Kp,Kv)

    dydt = zeros(2,1);

    dydt(1) = y(2);
    dydt(2) = ad(t) + Kp *( vd(t) - y(2) ) - Kv*y(2) + Kv*vd(t) - Kp*Kv*y(1) + Kv*Kp*pd(t);
    
end

