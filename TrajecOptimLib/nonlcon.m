function [c,ceq] = nonlcon(X,constants)

    if ~constants.uselogbar
        c = nonlcon_ineq(X, constants);
    else
        c = []; 
    end
    ceq = nonlcon_eq(X, constants);

end