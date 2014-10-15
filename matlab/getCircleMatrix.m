function circle = getCircleMatrix(x, y, r)

circle = eye(3,3);
circle(1,3) = -x;
circle(3,1) = -x;
circle(2,3) = -y;
circle(3,2) = -y;
circle(3,3) = -r^2 + x^2 + y^2;