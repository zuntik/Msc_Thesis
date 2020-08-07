function c = vehicle_collision_avoidance_bad(X1,X2,CONSTANTS)
    c = CONSTANTS.min_dist_intervehicles - min(sqrt(sum(BernsteinPow(X2-X1,2),2)));
end
