clear all
close all

%% Test the Basis equation

figure, hold on
fplot(@(times)arrayfun(@(t)BernsteinBasis(2,5,t),times),[0 1]);
fplot(@(x)(BernsteinEval([0,0,1,0,0,0],1,x)-0.01),[0 1])

%%

figure, hold on
p = [ 2 4 -1 2 5];

fplot(@(times)arrayfun(@(x) sum(arrayfun(@(i) (p(i)*BernsteinBasis(i-1,4,x)),1:5)), times),[0 1]);

BernsteinPlot(p-0.1,1,'PlotControlPoints',false);


%% Test for 2D
figure, hold on
p = [
0      0
cos(1) sin(1)
cos(2) sin(2)
cos(3) sin(3)
cos(4) sin(4)
cos(5) sin(5)
cos(6) sin(6)
0      0 
];

BernsteinPlot(p,1,'AddArrows',false);



%% Test for 3D
figure, hold on
view(3)
grid on
p = [
0 0      0
1 cos(1) sin(1)
2 cos(2) sin(2)
3 cos(3) sin(3)
4 cos(4) sin(4)
5 cos(5) sin(5)
6 cos(6) sin(6)
7 0      0 
];

BernsteinPlot(p,1,'PlotControlPoints',false);

%% Test toMonomial

p = [ 2 4 -1 2 5];
m = BernsteinToMon(p,1);

figure, hold on
fplot(@(x)polyval(m,x),[0 1]);
BernsteinPlot(p-0.1,1);

%% Test Deriv

p = [ 2 4 -1 2 5];
m = BernsteinToMon(p,1);

figure, hold on
fplot(@(x)polyval(polyder(m),x),[0 1]);
BernsteinPlot(BernsteinDeriv(p,1)-1,1);

%% Test AntiDeriv

p = [ 2 4 -1 2 5];
m = BernsteinToMon(p,2);

figure, hold on
fplot(@(x)polyval(polyint(m',2),x),[0 2]);
BernsteinPlot(BernsteinAntiDeriv(p,2,3)-.1,2);


%% Test Mul

p1 = [ 2 4 -1 2 5];
p2 = [ 5 1 8 ];
m1 = BernsteinToMon(p1,3);
m2 = BernsteinToMon(p2,3);

m = conv(m1,m2);
p = BernsteinMul(p1,p2);

figure, hold on
fplot(@(x)polyval(m,x),[0 3]);
BernsteinPlot(p-1,3);

%% Test Sum

p1 = [ 2 4 -1 2 5];
p2 = [ 5 1 8 ];
m1 = BernsteinToMon(p1,3);
m2 = BernsteinToMon(p2,3);

m = m1 + [ [ 0;  0];m2 ];
p = BernsteinSum(p1,p2);

figure, hold on
fplot(@(x)polyval(m,x),[0 3]);
BernsteinPlot(p-1,3);

%% Test Power

p = [ 2 4 -1 2 5];
m = BernsteinToMon(p,3);

p = BernsteinPow(p,2);
m = conv(m,m);

figure, hold on
fplot(@(x)polyval(m,x),[0 3]);
BernsteinPlot(p-1,3);

%% Test Integral

p = [ 2 4 -1 2 5];
m = BernsteinToMon(p,3);

disp(polyval(polyint(m'),3));
disp(BernsteinIntegr(p, 3));


%% Test Degree Elevation

p = [ 2 4 -1 2 5];
new_p = BernsteinDegrElev(p,10);
figure, hold on
BernsteinPlot(p-0.1, 5);
BernsteinPlot(new_p, 5);

%% Test Control Point evaluation
% objective is to test the curve in the times corresponding to the control
% points

p = [ 2 4 -1 2 5];
mat = BernsteinCtrlPntsEval(4);
figure, hold on
BernsteinPlot(p,1);
scatter(0:1/(length(p)-1):1, mat*p','filled');
