function constants = processconstants(constantsIn)
% function that adds more necessary information for the optimisation to run
% that includes calculating degree elevatiom matrices and setting default
% parameters

    constants = constantsIn;
    
    % isfied checks if a field in a structure is defined
    
    if isfield(constantsIn,'filled') && constantsIn.filled == true
        return
    end
    
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
        constants.min_dist_int_veh = 1;
    end
    
    if ~isfield(constants, 'obstacles')
        constants.obstacles = [];
    end
    
    if ~isfield(constants, 'obstacles_circles')
        constants.obstacles_circles = [];
    end
    
    if ~isfield(constants, 'plotboatsize')
        constants.plotboatsize = 1;
    end
    
    if ~isfield(constants, 'N')
        constants.N = 50;
    end
    
    if ~isfield(constants, 'useeqlogbar')
        constants.useeqlogbar = true;
    end
    
    if ~isfield(constants,'init_guess')
        constants.init_guess = @rand_init_guess;
    end
    
    constants.DiffMat = BernsteinDerivElevMat(constants.N,constants.T);
    constants.EvalMat = BernsteinCtrlPntsEval(constants.N);
    constants.BigElevMat = BernsteinDegrElevMat(constants.N,constants.N*10);
    constants.numvars = size(constants.xi,2);
    constants.Nv = size(constants.xi,1);%number of vehicles
        
    constants.filled = true;

end