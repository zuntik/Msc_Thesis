function xinit = rand_init_guess(constants)
    if ~isfield(constants,'statebounds') || isempty(constants.statebounds)
        xstatevars = rand(constants.numvars*(constants.N-1), 1);
    else
        statebounds = constants.statebounds;
        statebounds(isinf(statebounds)) = 1;
        xstatevars = (statebounds(2,:) - statebounds(1,:)) .* rand(constants.N-1, constants.numvars) + statebounds(1,:);
    end
    if ~isfield(constants,'inputbounds') || isempty(constants.inputbounds)
        xinputs = rand(constants.numinputs*(constants.N+1), 1);
    else
        xinputs = (constants.inputbounds(2,:) - constants.inputbounds(1,:)) .* rand(constants.N+1, constants.numinputs) + constants.inputbounds(1,:);
    end
    
    xinit = repmat([xstatevars(:);xinputs(:)],[1, constants.Nv]);
        
end
% function xinit = rand_init_guess(constants)
% %     xinit = rand((constants.numvars*(constants.N-1)+...
% %         constants.numinputs*(constants.N+1))*constants.Nv,1);
%     if ~isfield(constants,'statebounds') || isempty(constants.statebounds)
%         xstatevars = rand(constants.numvars*(constants.N-1)*constants.Nv,1);
%     else
%         xstatevars = (constants.statebounds(2,:) - constants.statebounds(1,:)) .* rand(constants.N-1, constants.numvars, constants.Nv) + constants.statebounds(1,:);
%     end
%     if ~isfield(constants,'inputbounds') || isempty(constants.inputbounds)
%         xinputs = rand(constants.numinputs*(constants.N+1)*constants.Nv,1);
%     else
%         xinputs = (constants.inputbounds(2,:) - constants.inputbounds(1,:)) .* rand(constants.N+1, constants.numinputs, constants.Nv) + constants.inputbounds(1,:);
%     end
%         
%     xinit = [ xstatevars(:); xinputs(:) ];
% end