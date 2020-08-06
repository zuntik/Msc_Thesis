Nc = 20;

circles = [ rand(Nc,1), rand(Nc,1), 0.1*rand(Nc,1)];

figure, axis equal, hold on

t = linspace(0,2*pi-0.01,1000);

indexes = circles_intersect(circles);
for i = 1:size(circles,1)
    c = [circles(i,3).*cos(t)+circles(i,1);circles(i,3).*sin(t)+circles(i,2)].';
    if ~indexes(i)
        plot(polyshape(c))
    else
        plot(c(:,1),c(:,2))
    end
end

