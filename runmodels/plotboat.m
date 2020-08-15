function plotboat(x,y,yaw,size)

    points = [
        size/2*cos(yaw + pi - pi/6 ), size/2*sin(yaw + pi - pi/6);
        size/2*cos(yaw + pi/6 ), size/2*sin(yaw + pi/6);
        size/1.5*cos(yaw), size/1.5*sin(yaw);
        size/2*cos(yaw - pi/6), size/2*sin(yaw-pi/6);
        size/2*cos(yaw - pi + pi/6),size/2*sin(yaw - pi + pi/6);
    ];

    points = points + [ x y];
    
    plot(polyshape(points));

end