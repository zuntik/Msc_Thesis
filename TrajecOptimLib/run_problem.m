function [xOut,J] = run_problem(constants)
% this function does each vehicle unconstrained then joins them all

    constants=processconstants(constants);

    xIn = reshape(constants.init_guess(constants), [], constants.Nv);

    options_fmincon = optimoptions(@fmincon,'Display','Iter','Algorithm','sqp','MaxFunctionEvaluations',Inf,'StepTolerance',eps,'MaxIterations',Inf);
    % options = optimoptions(@fmincon,'Display','Iter','Algorithm','sqp');
    options_fminsearch = optimset('Display','iter','Algorithm','sqp','MaxFunEvals',300000,'MaxIter',Inf);

    xOut_uc = zeros((constants.N-1)*constants.numvars+(constants.N+1)*constants.numinputs,constants.Nv);
    CONSTANTS2 = constants;
    CONSTANTS2.Nv=1;
    CONSTANTS2.obstacles = [];
    CONSTANTS2.obstacles_circles = [];
    [lb,ub] = getvariablebounds(CONSTANTS2);
    for i = 1:constants.Nv
        tic
        CONSTANTS2.xi = constants.xi(i,:);
        CONSTANTS2.xf = constants.xf(i,:);
        if constants.uselogbar
            if constants.useeqlogbar % use this if we want to put eq constraints on log bar
                [xOut_i,J] = fminsearch(@(x)costfun(x,CONSTANTS2),xIn(:,i),options_fminsearch);
%                 [xOut_i,J] = fmincon(@(x)costfun(x,CONSTANTS2),xIn(:,i),[],[]);
            else
                [xOut_i,J] = fmincon(@(x)costfun(x,CONSTANTS2),xIn(:,i),...
                    [],[],[],[],[],[],@(x)nonlcon(x,CONSTANTS2),options_fmincon);
            end
        else
            [xOut_i,J] = fmincon(@(x)costfun(x,CONSTANTS2),xIn(:,i),...
                [],[],[],[],lb,ub,@(x)nonlcon(x,CONSTANTS2),options_fmincon);
        end
        xOut_uc(:,i) = xOut_i;
    end
    elapsed_unconstrained = toc;
    disp(['time for unconstrained part: ', num2str(elapsed_unconstrained)])

    [lb,ub] = getvariablebounds(constants);
    if ~isempty(constants.obstacles) || constants.Nv>1 || ~isempty(constants.obstacles_circles)
        tic
        if constants.uselogbar
%             [xOut,J] = fminsearch(@(x)costfun(x,constants),xOut_uc,options_fminsearch);
%             [xOut,J] = fmincon(@(x)costfun(x,constants),xOut_uc(:),[],[]);
            if constants.useeqlogbar
                [xOut,J] = fminsearch(@(x)costfun(x,constants),xOut_uc(:),options_fminsearch);
%                 [xOut,J] = fmincon(@(x)costfun(x,constants),xOut_uc(:),[],[]);
            else
                [xOut,J] = fmincon(@(x)costfun(x,constants),xOut_uc(:),...
                    [],[],[],[],[],[],@(x)nonlcon(x,constants),options_fmincon);
            end
        else
            [xOut,J] = fmincon(@(x)costfun(x,constants),xOut_uc(:),[],[],[],[],lb,ub,@(x)nonlcon(x,constants),options_fmincon);
        end
        disp('Elapsed time for joined up problem:')
        elapsed_constrained = toc;
        disp(elapsed_constrained);
    else
        xOut = xOut_uc;     
    end
    disp(['time for unconstrained part: ', num2str(elapsed_unconstrained)])
%     disp(['time for constrained part: ', num2str(elapsed_constrained)])
    xOut = matrify(xOut,constants);

end