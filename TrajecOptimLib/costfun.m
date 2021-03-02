function J = costfun(X,constants)

    J = 0;
    if constants.uselogbar
        c = nonlcon_ineq(X, constants);
%         J = J + 10e-5*sum(logbarrierfunc(0.1, -c, constants.usesigma));
        J = J + 5*10e-2*sum(logbarrierfunc(0.4, -c, constants.usesigma));
        if constants.useeqlogbar
            ceq = nonlcon_eq(X, constants);
            %J=J+sum(logbarrierfunc(0.01,-abs(ceq),constants.usesigma));
            %J=J+sum(logbarrierfunc(0.01,1e-5-ceq.^2,constants.usesigma));
            %J=J+logbarrierfunc(0.1,-sum(1e-5-ceq.^2),constants.usesigma);
            %J = J + sum(ceq.^2);
            J=J+10e5*sum(...
                logbarrierfunc(0.1,-(ceq+10e-6),constants.usesigma)+...
                logbarrierfunc(0.1,(ceq+10e-6),constants.usesigma)...
                );
        end
    end
    
    X = matrify(X,constants);
    % dynamics
    J_dyn = zeros(constants.Nv,1);
    % calc cost for each vehicle and do sum
    for i = 1:constants.Nv
        J_dyn(i) = constants.costfun_single(X(:,:,i),constants);
    end    
    J = J + sum(J_dyn);

end