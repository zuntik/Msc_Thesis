function BernsteinPlot(p,T,varargin)
    hold on
    parser = inputParser;
    addParameter(parser,'PlotControlPoints',true);
    addParameter(parser,'AddArrows',false);
    parse(parser,varargin{:});
    
    [ n, dim ] = size(p);
    if n==1 && dim ~= 1
        p = p';
        dim = dim+n;
        n = dim - n;
        dim = dim- n;
    end
    
    if dim == 1
        fplot(@(times)arrayfun(@(t)(BernsteinEval(p,T,t)),times),[0 T]);
        if parser.Results.PlotControlPoints
            scatter(0:T/(n-1):T, p);
        end
    elseif dim == 2
        xt = @(times)arrayfun(@(t)BernsteinEval(p(:,1),T,t),times);
        yt = @(times)arrayfun(@(t)BernsteinEval(p(:,2),T,t),times);
        fplot(xt,yt,[0 T]);
        if parser.Results.PlotControlPoints
            scatter(p(:,1)',p(:,2)');
        end
        dp = BernsteinDeriv(p,T);
        for i = 0:9
            grad = BernsteinEval(dp,T,T*i/10);
            if parser.Results.AddArrows
                quiver(xt(T*i/10),yt(T*i/10),grad(1)/2,grad(2)/2,1,'linewidth',2,'color','r','MaxHeadSize',10);
            end
        end
        xlabel('Position x (m)');
        ylabel('Position y (m)');
    elseif dim == 3
        %xt = @(t) BernsteinEval(p(:,1),T,t);
        xt = @(times)arrayfun(@(t)BernsteinEval(p(:,1),T,t),times);
        %yt = @(t) BernsteinEval(p(:,2),T,t);
        yt = @(times)arrayfun(@(t)BernsteinEval(p(:,2),T,t),times);
        %zt = @(t) BernsteinEval(p(:,3),T,t);
        zt = @(times)arrayfun(@(t)BernsteinEval(p(:,3),T,t),times);
        fplot3(xt,yt,zt,[0 T]);
        view(3);
        if parser.Results.PlotControlPoints
            scatter3(p(:,1)',p(:,2)',p(:,3)');
        end
        view(3);
    end
end