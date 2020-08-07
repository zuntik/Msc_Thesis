function plot_xy(X,CONSTANTS)

    figure, axis equal, hold on
    for i = 1:CONSTANTS.Nv
        BernsteinPlot(X(:,1:2,i),CONSTANTS.T);
        [~,xy] = CONSTANTS.recoverxy(X(:,:,i),CONSTANTS.T);
        plot(xy(:,1),xy(:,2));
    end

    if ~isempty(CONSTANTS.obstacles) && true
        for i = 1:size(CONSTANTS.obstacles,3)
            plot(polyshape(CONSTANTS.obstacles(:,:,i)))
        end
    end

    if ~isempty(CONSTANTS.obstacles_circles)
        for i = 1:size(CONSTANTS.obstacles_circles,3)
            centrx = CONSTANTS.obstacles_circles(i,1);
            centry = CONSTANTS.obstacles_circles(i,2);
            r = CONSTANTS.obstacles_circles(i,3);
            plot(polyshape([(centrx+r*cos(0:0.01:2*pi));(centry+r*sin(0:0.01:2*pi))].'));
        end
    end
    
end