function constants = processconstants(constantsIn)

    constants = constantsIn;
    if isfield(constantsIn,'filled') && constantsIn.filled == true
        return
    end
    
    constants.filled = true;
    
    constants.DiffMat = BernsteinDerivElevMat(constants.N,constants.T);
    constants.EvalMat = BernsteinCtrlPntsEval(constants.N);
    constants.BigElevMat = BernsteinDegrElevMat(constants.N,constants.N*10);
    constants.numvars = size(constants.xi,2);
    constants.Nv = size(constants.xi,1);%number of vehicles
    
    if ~isfield(constants, 'uselogbar') 
        constants.uselogbar = false;
    end
    
    if ~isfield(constants, 'usesigma')
        constants.usesigma = true;
    end
    
    if ~isfield(constants, 'numinputs')
        constants.numinputs = 0;
    end
    
    if ~isfield(constants, 'min_dist_int_veh')
        constants.min_dist_int_veh = 0;
    end
    
    if ~isfield(constants, 'obstacles')
        constants.obstacles = [];
    end
    
    if ~isfield(constants, 'obstacles_circles')
        constants.obstacles_circles = [];
    end
    
end