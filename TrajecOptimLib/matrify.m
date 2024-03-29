function X = matrify(X,constants)
%     X = [constants.xi; reshape(X,[],constants.numvars); constants.xf];
    
%     X = cat(1,...
%         reshape(constants.xi.',[1,constants.numvars,constants.Nv]),...
%         reshape(X,[constants.N-1,constants.numvars,constants.Nv]),...
%         reshape(constants.xf.',[1,constants.numvars,constants.Nv]));

%     Xstatevars = X(1:constants.numvars*(constants.N-1)*constants.Nv);
%     Xinputs = X(constants.numvars*(constants.N-1)*constants.Nv+1:end);
%     
%     Xstatevars = cat(1,...
%         reshape(constants.xi.',[1,constants.numvars,constants.Nv]),...
%         reshape(Xstatevars,[constants.N-1,constants.numvars,constants.Nv]),...
%         reshape(constants.xf.',[1,constants.numvars,constants.Nv]));    
%     
%     Xinputs = reshape(Xinputs,[constants.N+1,constants.numinputs,constants.Nv]);
%     
%     X = cat(2, Xstatevars, Xinputs);
% 
    if ~iscolumn(X)
        return
    end
    constants = processconstants(constants);
    X = reshape(X, [], constants.Nv);
    Xstatevars = cat(1,...
        reshape(constants.xi.',[1,constants.numvars,constants.Nv]),...
        reshape(X(1:constants.numvars*(constants.N-1),:),[constants.N-1,constants.numvars,constants.Nv]),...
        reshape(constants.xf.',[1,constants.numvars,constants.Nv]));
    Xinputs = reshape(X(constants.numvars*(constants.N-1)+1:end,:),[constants.N+1,constants.numinputs,constants.Nv]);
    X = cat(2, Xstatevars, Xinputs);
    
end