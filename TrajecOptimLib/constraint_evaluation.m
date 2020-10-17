function constraint_evaluation(X,constants)
    
    constants = processconstants(constants);
    Xvars = X(2:end-1,1:constants.numvars,:);
    Xinputs = X(:,constants.numvars+1:end,:);
    X = [ Xvars(:);Xinputs(:)];
    [c,ceq] = nonlcon(X(:),constants);
    disp((c<0).')
    disp(norm(ceq))
    disp(ceq.')


end