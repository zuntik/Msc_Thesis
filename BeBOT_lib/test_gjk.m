shape1 = [ 0 0; 1 0; 1 1 ; 0 1 ];
shape2 =  shape1 + [ 0.5 0 ];
[dist, pt1, pt2, simplex] = gjk(shape1.', shape2.');
disp(['collision: ',num2str(simplex.collision),' distance: ',num2str(dist),' >0? ',num2str(dist>0)])
figure, hold on
plot(polyshape(shape1))
plot(polyshape(shape2))


shape3 = shape1;
shape4 = [ 0 0.5 ; 0.5 1; 1 0.5 ; 0.5 0];
[dist, pt1, pt2, simplex] = gjk(shape3.', shape4.');
disp(['collision: ',num2str(simplex.collision),' distance: ',num2str(dist),' >0? ',num2str(dist>0)])
figure, hold on
plot(polyshape(shape3))
plot(polyshape(shape4))

shape5 = shape1;
shape6 = shape4 + [ 0.5 0];
[dist, pt1, pt2, simplex] = gjk(shape5.', shape6.');
disp(['collision: ',num2str(simplex.collision),' distance: ',num2str(dist),' >0? ',num2str(dist>0)])
figure, hold on
plot(polyshape(shape5))
plot(polyshape(shape6))

shape7 = shape1;
shape8 = shape4 + [ 1 0];
[dist, pt1, pt2, simplex] = gjk(shape7.', shape8.');
disp(['collision: ',num2str(simplex.collision),' distance: ',num2str(dist),' >0? ',num2str(dist>0)])
figure, hold on
plot(polyshape(shape7))
plot(polyshape(shape8))