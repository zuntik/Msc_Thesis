function plot_xy(X,constants)
    
    constants=processconstants(constants);
    figure, axis equal, hold on, axis ij,camroll(90)
    for i = 1:constants.Nv
        BernsteinPlot(X(:,1:2,i),constants.T,'PlotControlPoints',false);
        [~,xy] = constants.recoverxy(X(:,:,i),constants);
        plot(xy(:,1),xy(:,2));
    end

    if ~isempty(constants.obstacles) && true
        for i = 1:size(constants.obstacles,3)
            plot(polyshape(constants.obstacles(:,:,i)))
        end
    end

    if ~isempty(constants.obstacles_circles)
        for i = 1:size(constants.obstacles_circles,3)
            centrx = constants.obstacles_circles(i,1);
            centry = constants.obstacles_circles(i,2);
            r = constants.obstacles_circles(i,3);
            plot(polyshape([(centrx+r*cos(0:0.01:2*pi));(centry+r*sin(0:0.01:2*pi))].'));
        end
    end
    
end