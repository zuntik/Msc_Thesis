clear all

P = rand(60,2);

tic
[min_val,indexes] = closest_points(P);
toc

figure, axis equal,hold on

scatter(P(:,1),P(:,2));
scatter(P(indexes,1),P(indexes,2),'filled');
plot(P(indexes,1),P(indexes,2))
