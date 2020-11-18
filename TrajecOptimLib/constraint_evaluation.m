function constraint_evaluation(X,constants)
    
    constants = processconstants(constants);
    Xvars = X(2:end-1,1:constants.numvars,:);
    Xinputs = X(:,constants.numvars+1:end,:);
    X = [ Xvars(:);Xinputs(:)];
    [c,ceq] = nonlcon(X(:),constants);
    disp((c<0).')
    disp(norm(ceq))
    disp(ceq.')

    if constants.Nv == 2
        figure
        X_diff = xOut(:,1:2,1)-xOut(:,1:2,2);
        X_diff_2 = BernsteinPow(X_diff,2);
        Mag = sqrt(sum(X_diff_2,2));
        BernsteinPlot(Mag,constants.T);
    end

    if false && (isfield(constants,'obstacles_circles') && ~isempty(constants.obstacles_circles))...
            && size(constants.obstacles_circles,3) == 1 && constants.Nv == 1
        figure, grid on
        BernsteinPlot(sum((xOut(:,1:2)-constants.obstacles_circles(1:2)).^2,2),constants.T);
    end

end