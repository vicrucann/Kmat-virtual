function v = homography_solve(S, UV)
n = size(UV, 2);
XY = zeros(size(UV));
for i = 1:n
    s = S{i};
    cs = getEllipseCenter(s);
    XY(:,i) = cs;
end
U = UV(1, :); V = UV(2,:); 
X = XY(1,:); Y = XY(2,:);
rows0 = zeros(3, n);
rowsXY = -[X; Y; ones(1,n)];
hx = [rowsXY; rows0; U.*X; U.*Y; U];
hy = [rows0; rowsXY; V.*X; V.*Y; V];
h = [hx hy];
if n == 4
    [U, ~, ~] = svd(h);
else
    [U, ~, ~] = svd(h, 'econ');
end
v = (reshape(U(:,9), 3, 3)).';
end