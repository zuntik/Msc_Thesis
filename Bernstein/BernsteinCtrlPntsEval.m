function [mat] = BernsteinCtrlPntsEval(order)
% the evaluation the the curve in the times corresponding to the control
% points are given by mat * cp, where the cp are displained in the
% following example form:
%     x1 y1 
%     x2 y2
%     x3 y3

    % define matrix - obtained from deCasteljau
    mat = (order:-1:0)'.^(order:-1:0) .* (0:order)'.^(0:order) ...
    .*  diag(rot90(pascal(order+1)))' ./ order^order;


end