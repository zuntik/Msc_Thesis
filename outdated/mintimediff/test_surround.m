a = 10e-8*rand(19,2);
a = [0 0; 1 0; 0.5 1000];


ca = convhull(a,'Simplify',true);
pca = polyshape(a(ca(1:end-1),:));

[centrx,centry] = centroid(pca);

r = max(norm(a(ca(1:end-1),:)-[centrx,centry]));
r = max(sqrt( sum((a(ca(1:end-1),:)-[centrx,centry]).^2,2)));

figure, hold on, axis equal
scatter(a(:,1),a(:,2));

plot(pca)

scatter(centrx,centry,'filled')
plot(polyshape([(centrx+r*cos(0:0.1:2*pi));(centry+r*sin(0:0.1:2*pi))].'))
