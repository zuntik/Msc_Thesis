T = 1;

p1 = [ 1 5 3 ];
p2 = [ 2 2 4 ];

figure, hold on
%BernsteinPlot(p1,T);
%BernsteinPlot(p2,T);

p1_2 = BernsteinPow(p1,2);
p2_2 = BernsteinPow(p2,2);

the_sqrt = @(t) sqrt(BernsteinEval(p1_2+p2_2,T,t));
fplot(the_sqrt,[0 T])

p1_2 = p1.^2;
p2_2 = p2.^2;

p3 = sqrt(p1_2+p2_2);
BernsteinPlot(p3,T);


% new example but simpler

clear p1 p2 p1_2 p2_2 p3


p1 = [ 1 5 ];
p2 = [ 2 2 ];

figure, hold on
%BernsteinPlot(p1,T);
%BernsteinPlot(p2,T);

p1_2 = BernsteinPow(p1,2);
p2_2 = BernsteinPow(p2,2);

the_sqrt = @(t) sqrt(BernsteinEval(p1_2+p2_2,T,t));
fplot(the_sqrt,[0 T])

p1_2 = p1.^2;
p2_2 = p2.^2;
p3 = sqrt(p1_2+p2_2);
BernsteinPlot(p3,T);


