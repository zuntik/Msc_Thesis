addpath('..\Bernstein\');
clear
close all

polygon = [ 0 0 ; 0 1; 1 1 ; 1 0];
cpts = [ 3.2 0; 2 1;  1 2; 0 3.2 ];

polygon = [ 0 0; 0 1; 1 1; 1 0 ]; % all in
cpts = [ 0.1 0.1; 0.1 0.9; 0.9 0.1; 0.9 0.9 ];

polygon = [ 0 0; 0 1; 1 1; 1 0 ]; % first and last in
cpts = [0.1 0.9; 0.1 1.9; 0.9 1.1; 0.9 0.1 ];

polygon = [ 0 0; 0 1; 1 1; 1 0 ]; %first and last are close
cpts = [ 0.5 1-10e-8; 0.5 1+10e-8 ];

polygon = [ 0 0 ; 0 1; 1 1 ; 1 0];
cpts = [ 0 2; 0.25 1.6; 0.5 1.6; 0.75 1.8; 1 2]; 

polygon = [ 0 0 ; 0 1; 1 1 ; 1 0];
cpts = [ 0 1.5; 0.25 1.1; 0.5 0.9; 0.75 0.5; 1 2];

polygon = [ 0 0 ; 0 1; 1 1 ; 1 0];
cpts = [ 0.1 0.9; 0.2 0.9; 0.5 1.1; 1.2 1.3 ];

polygon = [ -1 -1; -1 1; 1 1; 1 -1];
cpts=[
    0.6290    0.7082
   -2.3766    1.4319
   -2.1898   -1.0542
   -1.8519   -0.0107
   -0.2469    1.5922
    0.8617    0.4756
    1.7806    0.1821
   -0.0078   -0.8456
   -2.2561   -0.4415
   -0.9308    1.4700
   ];
cpts = 5*rand(10,2) - 2.5;
% cpts=[
%     0.6559   -1.5795
%    -0.7246    1.1289
%     2.4850   -0.6482
%    -1.3791    1.7078
%     0.7623    1.1711
%     0.5250    0.3551
%    -0.5638   -1.6157
%    -1.7891    2.2869
%    -2.3743   -1.1734
%    -0.3944    2.1229
% ];

tic
[dist, t, pt] = MinDistBernstein2Polygon(cpts.', polygon.');
dur1 = toc;
tic
[dist2,t2,pt2]= MinDistBernstein2Polygon_extended(cpts.', polygon.');
dur2 = toc;

figure, hold on
BernsteinPlot(cpts,1);

closest = BernsteinEval(cpts,1,t);
scatter(closest(1),closest(2))
scatter(pt(1),pt(2))
plot(polyshape(polygon))

figure, hold on
BernsteinPlot(cpts,1);
plot(polyshape(polygon))

if isempty(t2)
    disp('the curve is completely in the shape')
else
    closest = BernsteinEval(cpts,1,t2);
    scatter(pt2(1,:),pt2(2,:),'filled')
end

disp(['dur1: ',num2str(dur1),' dur2: ',num2str(dur2)]);
