clear 
addpath('Bernstein')

n=100;
a = 100*rand(n+1,2);

% tic
% b = BernsteinDegrElev(a,10000);
% c = BernsteinEval(a,1,linspace(0,1,10001));
% toc

tic
elevmat = BernsteinDegrElevMat(n,10000);
evalmat = BernsteinEvalMat(n,1,linspace(0,1,10001));
toc

tic
b = elevmat*a;
c = evalmat*a;
toc
%%
figure
BernsteinPlot(a,1);

scatter(b(:,1),b(:,2),'filled')
scatter(c(:,1),c(:,2),'filled')
