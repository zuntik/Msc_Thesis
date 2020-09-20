function [xOut,J] = run_problem(CONSTANTS)
% this function does each vehicle unconstrained then joins them all

    
    xIn = CONSTANTS.init_guess(CONSTANTS);

    options_fmincon = optimoptions(@fmincon,'Display','Iter','Algorithm','sqp','MaxFunctionEvaluations',300000,'StepTolerance',eps,'MaxIterations',Inf);
    %options_fminsearch = optimoptions(@fminsearch,'Display','iter','Algorithm','sqp','MaxFunEvals',300000,'MaxIter',Inf);
    % options = optimoptions(@fmincon,'Display','Iter','Algorithm','sqp');

    %xOut_uc = zeros(CONSTANTS.N-1,CONSTANTS.numvars,CONSTANTS.Nv);
    xOut_uc = zeros((CONSTANTS.N-1)*CONSTANTS.numvars+(CONSTANTS.N+1)*CONSTANTS.numinputs,CONSTANTS.Nv);
    CONSTANTS2 = CONSTANTS;
    CONSTANTS2.Nv=1;
    CONSTANTS2.obstacles = [];
    CONSTANTS2.obstacles_circles = [];
    for i = 1:CONSTANTS.Nv
        tic
        CONSTANTS2.xi = CONSTANTS.xi(i,:);
        CONSTANTS2.xf = CONSTANTS.xf(i,:);
        if CONSTANTS.uselogbar
            %[xOut_i,J] = fminsearch(@(x)costfun(x,CONSTANTS2),xIn);
            [xOut_i,J] = fmincon(@(x)costfun(x,CONSTANTS2),xIn(:,i),[],[]);
        else
            [xOut_i,J] = fmincon(@(x)costfun(x,CONSTANTS2),xIn(:,i),...
                [],[],[],[],[],[],@(x)nonlcon(x,CONSTANTS2),options_fmincon);
        end
        xOut_uc(:,i) = xOut_i;
        %xOut_uc(:,:,i) = reshape(xOut_i,(CONSTANTS.N-1),CONSTANTS.numvars);
        disp(['time for vehicle ', num2str(i), ':'])
        toc
    end
    
    xOut_uc_states = xOut_uc(1:CONSTANTS.numvars*(CONSTANTS.N-1),:);
    xOut_uc_inputs= xOut_uc(CONSTANTS.numvars*(CONSTANTS.N-1)+1:end,:);
    xOut_uc = [xOut_uc_states(:); xOut_uc_inputs];
    
    if ~isempty(CONSTANTS.obstacles) || CONSTANTS.Nv>1 || ~isempty(CONSTANTS.obstacles_circles)
        tic
        if CONSTANTS.uselogbar
            %[xOut,J] = fminsearch(@(x)costfun(x,CONSTANTS),xOut_uc,options_fminsearch);
            [xOut,J] = fmincon(@(x)costfun(x,CONSTANTS),xOut_uc,[],[]);
        else
            [xOut,J] = fmincon(@(x)costfun(x,CONSTANTS),xOut_uc,[],[],[],[],[],[],@(x)nonlcon(x,CONSTANTS),options_fmincon);
        end
        disp('Elapsed time for joined up problem:')
        toc
    else
        xOut = xOut_uc;     
    end

    xOut = matrify(xOut,CONSTANTS);

    
end