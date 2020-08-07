function X = matrify(X,CONSTANTS)
%     X = [CONSTANTS.xi; reshape(X,[],CONSTANTS.numvars); CONSTANTS.xf];
    
%     X = cat(1,...
%         reshape(CONSTANTS.xi.',[1,CONSTANTS.numvars,CONSTANTS.Nv]),...
%         reshape(X,[CONSTANTS.N-1,CONSTANTS.numvars,CONSTANTS.Nv]),...
%         reshape(CONSTANTS.xf.',[1,CONSTANTS.numvars,CONSTANTS.Nv]));

    Xstatevars = X(1:CONSTANTS.numvars*(CONSTANTS.N-1)*CONSTANTS.Nv);
    Xinputs = X(CONSTANTS.numvars*(CONSTANTS.N-1)*CONSTANTS.Nv+1:end);
    
    Xstatevars = cat(1,...
        reshape(CONSTANTS.xi.',[1,CONSTANTS.numvars,CONSTANTS.Nv]),...
        reshape(Xstatevars,[CONSTANTS.N-1,CONSTANTS.numvars,CONSTANTS.Nv]),...
        reshape(CONSTANTS.xf.',[1,CONSTANTS.numvars,CONSTANTS.Nv]));    
    
    Xinputs = reshape(Xinputs,[CONSTANTS.N+1,CONSTANTS.numinputs,CONSTANTS.Nv]);
    
    X = cat(2, Xstatevars, Xinputs);
    
end