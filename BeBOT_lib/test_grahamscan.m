clear, close all

A = [ 0 0 ; -1 0; 1 0 ; 0 1; 1 1 ];
A = [ 1:10;1:10].';
A = [1:10; zeros(1,10)].';
A = [zeros(1,10);1:10].';
A = [ 10:-1:1;10:-1:1 ].';
A= [rand(100,1),rand(100,1)];
A1 = [ 0 0; 0 1; 1 1; 1 0 ];
A2 = [ 0.1 0.1; 0.1 0.9; 0.9 0.9; 0.9 0.1 ];
A2 = [ 0.1 0.1; 0.1 0.9; 0.9 0.1; 0.9 0.9 ];
A = [ A1; A2];
A = [A2;A1];

indexes = grahamscan(A);

convexhull=A(indexes,:);
figure
hold on
scatter(A(:,1),A(:,2));
scatter(convexhull(:,1),convexhull(:,2),'filled')