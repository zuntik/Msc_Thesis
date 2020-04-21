function J = costfun_tt(x,T,Kp,Kv)
% cost function for a 1-D vehicle with cost given by the integral of v^2 
% INPUT
% x: control points
% T: time domain
% OUTPUT
% J: cost

y0 = [
    x(1) % initial position
    (x(2) - x(1)) * (length(x)-1)/T % initial speed
    0 % initial integral of v^2
];

dx  = BernsteinDeriv(x, T);
ddx = BernsteinDeriv(dx,T);

pd = @(t) BernsteinEval(x,  T,t);
vd = @(t) BernsteinEval(dx, T,t);
ad = @(t) BernsteinEval(ddx,T,t);

[~,y] = ode45(@(t,y) TrajectoryTrackingODE(t,y,pd,vd,ad,Kp,Kv), [0 T], y0);

J = y(end,3);

end