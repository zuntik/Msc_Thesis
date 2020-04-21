function dydt = TrajectoryTrackingODE(t,y,pd,vd,ad,Kp,Kv)

    dydt = zeros(3,1);

    dydt(1) = y(2);
    dydt(2) = ad(t) + Kp *( vd(t) - y(2) ) - Kv*y(2) + Kv*vd(t) - Kp*Kv*y(1) + Kv*Kp*pd(t);
    dydt(3) = y(2)^2;

end