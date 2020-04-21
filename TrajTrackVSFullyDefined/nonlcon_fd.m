function [c, ceq] = nonlcon_fd(x,n,T, b, k)

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
p = x(1:n+1);
v = x(n+2: 2*n+2);
u = x(2*n+3: 3*n+3);

% calculate derivatives

dv = BernsteinDeriv(v,T);

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