function [d,coor_min] = CurveDist2p(p,T,p0)
    % assumes that pol is 2 dimensional
    
    % center pol to p0 as origin
    p0 = p0(:);
    shift_p = BernsteinSum(p,-p0');
    
    ptimesdp = sum(BernsteinMul(BernsteinDeriv(shift_p,T),shift_p),2);
    
    min_t = roots(BernsteinToMon(ptimesdp,T)');
    real_min_t = min_t(min_t == real(min_t) & min_t>0 & min_t<T );
    
    real_min_t = [ 0 , real_min_t', T];
    
    coor_roots = BernsteinEval(shift_p,T,real_min_t);
    
    [d,idx] = min(arrayfun(@(i)norm(coor_roots(:,i)), 1:size(coor_roots,2)));
    
    if nargout > 1
        coor_min = BernsteinEval(p,T,real_min_t(idx));
    end
    
end

