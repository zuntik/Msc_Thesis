function [X,J] = run_problem_progressive_n(constants)
% this function does each vehicle unconstrained then joins them all

    constants=processconstants(constants);

    X = constants.init_guess(constants);
    
    options_fmincon = optimoptions(@fmincon,'Display','Iter','Algorithm','sqp','MaxFunctionEvaluations',Inf,'StepTolerance',eps,'MaxIterations',Inf);
    %options = optimoptions(@fmincon,'Display','Iter','Algorithm','sqp');
    options_fminsearch = optimset('Display','iter','Algorithm','sqp','MaxFunEvals',300000,'MaxIter',Inf);

    [lb,ub] = getvariablebounds(constants);
    
    J_old = Inf;
    while true
        tic
        if constants.uselogbar
            if constants.useeqlogbar
                [X,J] = fminsearch(@(x)costfun(x,constants),X,options_fminsearch);
                %[X,J] = fmincon(@(x)costfun(x,constants),X,[],[]);
            else
                [X,J] = fmincon(@(x)costfun(x,constants),X,...
                    [],[],[],[],[],[],@(x)nonlcon(x,constants),options_fmincon);
            end
        else
            [X,J] = fmincon(@(x)costfun(x,constants),X,[],[],[],[],lb,ub,@(x)nonlcon(x,constants),options_fmincon);
        end
        toc
        if abs((J_old-J)/J)<0.01
            break
        end
        
        J_old = J;
        X = increaseorder(X, constants, constants.N+10);
        constants.N = constants.N+10;
        [lb,ub] = getvariablebounds(constants);
        constants.filled = false;
        constants = processconstants(constants);
        
    end

    X = matrify(X,constants);

end