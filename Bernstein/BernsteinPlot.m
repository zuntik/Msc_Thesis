function BernsteinPlot(p,T)
    [ n, dim ] = size(p);
    if dim == 1
        fplot(@(t)(BernsteinEval(p,T,t)),[0 T]);
        scatter(0:T/(n-1):T, p);
    elseif dim == 2
        xt = @(t) BernsteinEval(p(:,1),T,t);
        yt = @(t) BernsteinEval(p(:,2),T,t);
        fplot(xt,yt,[0 T]);
        scatter(p(:,1)',p(:,2)');
        dp = BernsteinDeriv(p,T);
        
        for i = 0:9
            grad = BernsteinEval(dp,T,T*i/10);
            quiver(xt(T*i/10),yt(T*i/10),grad(1)/2,grad(2)/2,1,'linewidth',2,'color','r','MaxHeadSize',10);
        end
        xlabel('Position x (m)');
        ylabel('Position y (m)');
    end
end