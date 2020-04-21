function [c, ceq] = noncolon_unicycle(x,n,T)


persistent old_n old_T  ctrl_points_mat degr_elev_mat_v degr_elev_mat_dv degr_elev_mat_u

if  isempty(old_n) || isempty(old_T) || old_n ~= n || old_T ~= T
    old_n = n;
    old_T = T;
    ctrl_points_mat = BernsteinCtrlPntsEval(2*n);
    degr_elev_mat_v = BernsteinDegrElevMat(n,2*n);
    degr_elev_mat_dv = BernsteinDegrElevMat(n-1,2*n);
    degr_elev_mat_u = BernsteinDegrElevMat(n,2*n);
end

% extract the variables
px = x(0*(n+1)+1:1*(n+1));
py = x(1*(n+1)+1:2*(n+1));
v  = x(2*(n+1)+1:3*(n+1));
r  = x(3*(n+1)+1:4*(n+1));

% calculate derivatives
dx = BernsteinDeriv(px,T);
dy = BernsteinDeriv(py,T);
ddx = BernsteinDeriv(dx,T);
ddy = BernteinsDeriv(dy,T);

%% calculate the squares
v_2 = BernsteinPow(v,2);
dx_2 = BernsteinPow(dx,2);
dy_2 = BernsteinPow(dy,2);

%% conditions - match control points on left and right
cond1 = v_2 - dx_2 - dy_2;

right_side = BernsteinSum(BernsteinMul(dx,ddy), -BernsteinMul(ddx,dy));
left_side = BernsteinMul(v_2,r);
cond2 = left_side-right_side;

% degree elevation
dv_elev = degr_elev_mat_dv * dv;
v_elev = degr_elev_mat_v * v;
u_elev = degr_elev_mat_u * u;

% Evalutaion of state variables and derivatives in control points
v_dot_hat = ctrl_points_mat * dv_elev;
v_hat = ctrl_points_mat * v_elev;
u_hat = ctrl_points_mat * u_elev;


% equate on the dynamics equation

condition = v_dot_hat + b .* v_hat .* abs(v_hat) + k.* u_hat;

ceq = norm(condition);
c = [];


end
