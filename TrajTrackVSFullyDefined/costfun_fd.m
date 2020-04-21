function J = costfun_fd(x,n,T)
% calculates cost by performing summation of control points

% extract the variables
p = x(1:n+1);
v = x(n+2: 2*n+2);
u = x(2*n+3: 3*n+3);

v2 = BernsteinPow(v,2);

J = BernsteinIntegr(v2,T);
% integration is simply the sum of the control points


end