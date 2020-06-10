clear
close  all
% this will use my gtk2 alg and Calvin's too
addpath('Bernstein')
addpath('BeBOT_lib')

obstacle = [ 2 2; 3 1; 4 1.5; 3.5 3; 2.5 2.7];


cpts = [ 2 1; 3 1.3; 3.3 2.5; 3.4 2.6; 3.9 2; 4.5 3];
N = size(cpts,1)-1;

figure,hold on
plot(polyshape(obstacle))
BernsteinPlot(cpts, 1);

[dist, t, pt] = MinDistBernstein2Polygon_extended(cpts.', obstacle.');

if dist<0
    scatter(pt(1,:),pt(2,:),'filled')
    if length(t) == 2
        cptsleftcut = deCasteljau(cpts.',t(1));
        cptsrightcut = deCasteljau(cptsleftcut(:,N+1:end),t(2));
        cpts_inside = cptsrightcut(:,1:N+1);
    else
        cpts_inside = pt;
    end
    p2 = ConvexPolygon(cpts_inside.');
    BernsteinPlot(cpts_inside.',1);
    p1 = ConvexPolygon(obstacle);    
    [~,simplex] = gjk2(p1,p2);
    support=@(d) p1.support(d)-p2.support(-d);
    vec_up = directed_epa(support, simplex, [0;1]);
    vec_side = directed_epa(support, simplex, [1;0]);
    
    disp(vec_up)
    disp(vec_side)
    
end

plot(polyshape(p1.matrix))
plot(polyshape(cpts_inside.'))
plot(polyshape(cpts_inside.'+vec_up.'))
plot(polyshape(cpts_inside.'+vec_side.'))

for i = 1:3
    p(i,:) = simplex{i}.P.'; %#ok<SAGROW>
end
plot(polyshape(p))
scatter(0,0,'filled')