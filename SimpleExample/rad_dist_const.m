function [c, ceq] = rad_dist_const(p,T,centre,min_dist)
    ceq = 0;
    for i = 1:size(centre,1)
        d(i) = CurveDist2p(p,T,centre(i,:));
        c(i) = min_dist(i) - d(i);
    end
    c = max(c);
end