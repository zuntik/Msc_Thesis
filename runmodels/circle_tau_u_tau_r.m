
fsolve(@func_for_circle, zeros(1, 4))
disp('u, v, tau_u, tau_r')

function F = func_for_circle(x)
% objective of fsolve is to make F = 0

T = 70;
radius = 10;
r = 2*pi/T;

mass = 17.0;
% I_z = 1;
X_dot_u = -20;
Y_dot_v = -30;%-1.3175;
% N_dot_r = -8.69;% -0.5;
X_u = -0.2;
Y_v = -50;
N_r = -4.14; %-0.1
X_uu = -25; 
Y_vv = -0.01;%-101.2776;
N_rr = -6.23; %-21

%%%%%% masses
m_u = mass - X_dot_u;
m_v = mass - Y_dot_v;
% m_r = I_z - N_dot_r;
m_uv = m_u - m_v;

u = x(1);
v = x(2);
tau_u = x(3);
tau_r = x(4);

d_u = -X_u -X_uu*abs(u);
d_v = -Y_v -Y_vv*abs(v);
d_r = -N_r -N_rr*abs(r);

F(1) = tau_u + m_v * v * r - d_u * u;
F(2) = m_u*u*r + d_v * v;
F(3) = tau_r + m_uv * u*v - d_r * r;
F(4) = u^2 + v^2 - radius^2*r^2;

end
