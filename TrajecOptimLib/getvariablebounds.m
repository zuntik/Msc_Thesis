function [lb, ub] = getvariablebounds(constants)
    
    if ~isfield(constants,'statebounds') || isempty(constants.statebounds)
        lb_s = Inf*ones((constants.N-1)*constants.numvars,1);
        ub_s = Inf*ones((constants.N-1)*constants.numvars,1);
    else
        lb_s = repmat(constants.statebounds(1,:), constants.N-1, 1);
        ub_s = repmat(constants.statebounds(2,:), constants.N-1, 1);
    end
    
%     lb_s = repmat(lb_s(:), [constants.Nv 1]);
%     ub_s = repmat(ub_s(:), [constants.Nv 1]);
    
    if ~isfield(constants,'inputbounds') || ...
        isempty(constants.inputbounds) || ...
        constants.numinputs == 0      
        lb_i = Inf*ones((constants.N+1)*constants.numinputs,1);
        ub_i = Inf*ones((constants.N+1)*constants.numinputs,1);
    else
        lb_i = repmat(constants.inputbounds(1,:), constants.N+1, 1);
        ub_i = repmat(constants.inputbounds(2,:), constants.N+1, 1);        
    end
    
    lb = repmat([lb_s(:);lb_i(:)], [constants.Nv,1]);
    ub = repmat([ub_s(:);ub_i(:)], [constants.Nv,1]);
    
    if (~isfield(constants,'statebounds')||isempty(constants.statebounds)) && ...
       (~isfield(constants,'inputbounds')||isempty(constants.inputbounds))
        lb = [];
        ub = [];
    end
    
end
% function [lb, ub] = getvariablebounds(constants)
%     
%     if ~isfield(constants,'statebounds') || isempty(constants.statebounds)
%         lb_s = Inf*ones((constants.N-1)*constants.numvars,1);
%         ub_s = Inf*ones((constants.N-1)*constants.numvars,1);
%     else
%         lb_s = repmat(constants.statebounds(1,:), constants.N-1, 1);
%         ub_s = repmat(constants.statebounds(2,:), constants.N-1, 1);
%     end
%     
%     lb_s = repmat(lb_s(:), [constants.Nv 1]);
%     ub_s = repmat(ub_s(:), [constants.Nv 1]);
%     
%     if ~isfield(constants,'inputbounds') || ...
%         isempty(constants.inputbounds) || ...
%         constants.numinputs == 0      
%         lb_i = Inf*ones((constants.N+1)*constants.numinputs,1);
%         ub_i = Inf*ones((constants.N+1)*constants.numinputs,1);
%     else
%         lb_i = repmat(constants.inputbounds(1,:), constants.N+1, 1);
%         ub_i = repmat(constants.inputbounds(2,:), constants.N+1, 1);        
%     end
%     
%     lb_i = repmat(lb_i(:), [constants.Nv 1]);
%     ub_i = repmat(ub_i(:), [constants.Nv 1]);
%     
%     lb = [lb_s; lb_i];
%     ub = [ub_s; ub_i];
%     
%     if (~isfield(constants,'statebounds')||isempty(constants.statebounds)) && ...
%        (~isfield(constants,'inputbounds')||isempty(constants.inputbounds))
%         lb = [];
%         ub = [];
%     end
%     
% end