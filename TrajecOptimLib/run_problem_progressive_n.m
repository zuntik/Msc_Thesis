function [X,J, X_history, J_history, Times_history] = run_problem_progressive_n(constants)
% this function does each vehicle unconstrained then joins them all

    constants=processconstants(constants);

    X = constants.init_guess(constants);

    options_fmincon = optimoptions(@fmincon,'Display','Iter','Algorithm','sqp','MaxFunctionEvaluations',Inf,'StepTolerance',eps,'MaxIterations',Inf);
    %options = optimoptions(@fmincon,'Display','Iter','Algorithm','sqp');
    options_fminsearch = optimset('Display','iter','Algorithm','sqp','MaxFunEvals',300000,'MaxIter',Inf);

    [lb,ub] = getvariablebounds(constants);

    J_history = zeros(1,8);
    Times_history = zeros(1,size(J_history,2));
    X_history = cell(1,size(J_history,2));
    
    for currun = 1:length(J_history)
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
        Times_history(currun) = toc;
        J_history(currun) = J;
        X_history{currun} = X;
        if currun~=1 && abs((J_history(currun-1)-J)/J_history(currun-1))<0.01
            break
        end
        
        figure
        plot_xy(X, constants);
        X = increaseorder(X, constants, constants.N+10);
        constants.N = constants.N+constants.Nincrement;
        [lb,ub] = getvariablebounds(constants);
        constants.filled = false;
        constants = processconstants(constants);

    end
    disp('times for each run')
    disp(Times_history);
    disp('costs');
    disp(J_history);
    X = matrify(X,constants);
    figure
    
end