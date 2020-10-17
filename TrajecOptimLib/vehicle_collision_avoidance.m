function c = vehicle_collision_avoidance(X1,X2,constants)
    c = constants.min_dist_int_veh - MinDistBernstein2Polygon((X2-X1).',[0;0]);
end