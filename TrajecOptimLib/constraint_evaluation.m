function constraint_evaluation(X,CONSTANTS)

    Xvars = X(2:end-1,1:CONSTANTS.numvars,:);
    Xinputs = X(:,CONSTANTS.numvars+1:end,:);
    X = [ Xvars(:);Xinputs(:)];
    [c,ceq] = nonlcon(X(:),CONSTANTS);
    disp((c<0).')
    disp(norm(ceq))
    disp(ceq.')


end