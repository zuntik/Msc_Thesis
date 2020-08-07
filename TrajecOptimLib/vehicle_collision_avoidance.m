function c = vehicle_collision_avoidance(X1,X2,CONSTANTS)
    c = CONSTANTS.min_dist_intervehicles - MinDistBernstein2Polygon((X2-X1).',[0;0]);
end