function center = getEllipseCenter(ellipse)

%center = -inv(ellipse(1:2, 1:2))*ellipse(1:2,3);

A = ellipse(1,1);
B = ellipse(1,2);
C = ellipse(2,2);
D = ellipse(1,3);
E = ellipse(2,3);

center = [(B*E - C*D) / (A*C - B^2); (B*D - A*E) / (A*C - B^2)];