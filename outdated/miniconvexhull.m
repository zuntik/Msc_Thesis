close all, clear all
% cpts = rand(10,2)*5;
% cpts =[
%     0.4669    2.5271
%     1.5368    3.8071
%     2.2803    3.1553
%     0.5083    0.4495
%     4.9769    0.4043
%     1.6605    3.8862
%     1.4867    4.5257
%     0.3102    2.6689
%     1.4912    0.5458
%     0.2318    4.1290
% ];
cpts =[
    0.4669    2.5271
    1.4298    3.6791
    2.1316    3.2857
    1.0399    1.2612
    3.1895    0.4224
    3.3187    2.1452
    1.5910    4.1420
    1.1338    3.9687
    0.5464    2.2443
    1.3653    0.9041
    0.2318    4.1290
];
cpts = [ 1 1; 1 0; 0 0.8; 10 0.5];
%cpts = [ 1 4; 1 3; 1 2; 1 1; 1 0; 0 0.5; 10 0.5];
%cpts = BernsteinDegrElev(cpts,7);
%cpts = [ 2 1; 3 1.3; 3.3 2.5; 3.4 2.6; 3.9 2; 4.5 3]+ [-2.5 -2];
N = size(cpts,1)-1;

% syms t
% B = bernsteinMatrix(N, t);
% bezierCurve = simplify(B*cpts);

figure,hold on%, axis equal
%fplot(bezierCurve(1), bezierCurve(2), [0, 1])

% curvepoints = bernsteinMatrix(N,linspace(0,1,N+1))*cpts;
%cuvepoints = BernsteinPolyStable(cpts.',linspace(0,1,N+1)).';
curvepoints = BernsteinCtrlPntsEval(N)*cpts;


%BernsteinPlot(cpts, 1);
manypoints = bernsteinMatrix(N,linspace(0,1,100000))*cpts;
% manypoints = BernsteinPolyStable(cpts.',linspace(0,1,100000)).';
plot(manypoints(:,1),manypoints(:,2));
scatter(cpts(:,1),cpts(:,2));
scatter(curvepoints(:,1),curvepoints(:,2),'filled');

for i = 1:N
    vertices = [cpts(i:i+1,:);curvepoints(i:i+1,:)];
    vertices=vertices(convhull(vertices,'simplify',true),:);
    plot(polyshape(vertices(1:end-1,:)));
    %scatter(vertices(:,1),vertices(:,2))
end

% dX = BernsteinDeriv(cpts,1);
% figure
% BernsteinPlot(cpts(:,2),1)
% scatter(linspace(0,1,N+1),curvepoints(:,2),'filled')
% plot([0 0.1],cpts(1:2,2));
% plot([0.9 1],cpts(end-1:end,2));
% 
% figure
% BernsteinPlot(cpts(:,1),1)
% scatter(linspace(0,1,N+1),curvepoints(:,1),'filled')
% plot([0 0.1],cpts(1:2,1));
% plot([0.9 1],cpts(end-1:end,1));
% 
% figure,hold on
% disp(cpts(1:5,:))
% BernsteinPlot(cpts(1:5,:),1)
% plot(cpts(1:2,1),cpts(1:2,2))
% p = BernsteinEval(cpts(1:5,:),1,0.2);
% scatter(p(1),p(2),'filled')
% plot([p(1) cpts(2,1)],[p(2) cpts(2,2)])
% plot([p(1) cpts(1,1)],[p(2) cpts(1,2)])