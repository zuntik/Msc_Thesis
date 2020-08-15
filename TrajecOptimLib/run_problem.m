function xOut = run_problem(CONSTANTS)
% this function does each vehicle unconstrained then joins them all
    
    tic
    xIn = CONSTANTS.init_guess(CONSTANTS);

    options = optimoptions(@fmincon,'Display','Iter','Algorithm','sqp','MaxFunctionEvaluations',300000,'StepTolerance',eps,'MaxIterations',Inf);
    % options = optimoptions(@fmincon,'Display','Iter','Algorithm','sqp');

    %xOut_uc = zeros(CONSTANTS.N-1,CONSTANTS.numvars,CONSTANTS.Nv);
    xOut_uc = zeros((CONSTANTS.N-1)*CONSTANTS.numvars+(CONSTANTS.N+1)*CONSTANTS.numinputs,CONSTANTS.Nv);
    CONSTANTS2 = CONSTANTS;
    CONSTANTS2.Nv=1;
    CONSTANTS2.obstacles = [];
    CONSTANTS2.obstacles_circles = [];
    for i = 1:CONSTANTS.Nv
        CONSTANTS2.xi = CONSTANTS.xi(i,:);
        CONSTANTS2.xf = CONSTANTS.xf(i,:);
        xOut_i = fmincon(@(x)costfun(x,CONSTANTS2),xIn,...
            [],[],[],[],[],[],@(x)nonlcon(x,CONSTANTS2),options);
        xOut_uc(:,i) = xOut_i;
        %xOut_uc(:,:,i) = reshape(xOut_i,(CONSTANTS.N-1),CONSTANTS.numvars);
    end
    
    xOut_uc_states = xOut_uc(1:CONSTANTS.numvars*(CONSTANTS.N-1),:);
    xOut_uc_inputs= xOut_uc(CONSTANTS.numvars*(CONSTANTS.N-1)+1:end,:);
    xOut_uc = [xOut_uc_states(:); xOut_uc_inputs];
    
    if ~isempty(CONSTANTS.obstacles) || CONSTANTS.Nv>1 || ~isempty(CONSTANTS.obstacles_circles)
        xOut = fmincon(@(x)costfun(x,CONSTANTS),xOut_uc,[],[],[],[],[],[],@(x)nonlcon(x,CONSTANTS),options);
    else
        xOut = xOut_uc;
    end

    xOut = matrify(xOut,CONSTANTS);
    
    toc
end