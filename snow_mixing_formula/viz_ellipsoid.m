clear

a = 1;
b = a;
c = 5;
[X,Y,Z] = ellipsoid(0,0,0,a,b,c);
figure(Position=[5.8966e+03 130.2000 560 420])
surf(X,Y,Z);
axis equal