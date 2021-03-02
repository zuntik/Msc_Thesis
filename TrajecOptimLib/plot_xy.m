function plot_xy(X,constants)

    constants=processconstants(constants);
    X=matrify(X, constants);

%     figure, axis equal, hold on, axis ij,camroll(90)
    axis equal
    for i = 1:constants.Nv
        BernsteinPlot(X(:,[2,1],i),constants.T,'PlotControlPoints',false);
        [~,~, tenpoints] = constants.recoverxy(X(:,:,i),constants);
%         [~,xy, tenpoints] = constants.recoverxy(X(:,:,i),constants);
%         plot(xy(:,2),xy(:,1));
        for v = 1:size(constants.xi, 1)
            for b = 1:10
                plotboat(tenpoints(b,2),tenpoints(b,1),pi/2-tenpoints(b,3),constants.plotboatsize);
            end
        end
    end

    if ~isempty(constants.obstacles) && true
        for i = 1:size(constants.obstacles,3)
            plot(polyshape(constants.obstacles(:,[2,1],i)))
        end
    end

    if ~isempty(constants.obstacles_circles)
        for i = 1:size(constants.obstacles_circles,3)
            centrx = constants.obstacles_circles(i,1);
            centry = constants.obstacles_circles(i,2);
            r = constants.obstacles_circles(i,3);
            plot(polyshape([(centry+r*cos(0:0.01:2*pi));(centrx+r*sin(0:0.01:2*pi))].'));
        end
    end

end