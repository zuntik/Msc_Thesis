function c = vehicle_collision_avoidance_isaac(X1,X2,CONSTANTS)
    c = CONSTANTS.min_dist_intervehicles -  min(sqrt(sum((CONSTANTS.BigElevMat*(X1-X2)).^2,2)));
end
