function constraint_evaluation(X,CONSTANTS)

    X = X(2:end-1,:,:);
    [c,ceq] = nonlcon(X(:),CONSTANTS);
    disp((c<0).')
    disp(norm(ceq))
    disp(ceq.')


end