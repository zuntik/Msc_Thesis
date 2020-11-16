function c = vehicle_collision_avoidance_fast(X1,X2,constants)
    c = constants.min_dist_int_veh -  min(sqrt(sum((X1-X2).^2,2)));
end
