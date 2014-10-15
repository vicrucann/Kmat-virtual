function RMS = calc_reprojection(UV, H, S)
n = size(S,3);
err = zeros(1,n);
h_inv = inv(H);
% % d(UV, HCS)
for i=1:n
    s = S{i};
    cs = getEllipseCenter(s);
    hs = h_inv'*s*h_inv;
    chs = getEllipseCenter(hs);
    hcs = reshape(H, 3, 3)*[cs(1) cs(2) 1]';
    hcs = hcs / hcs(3);
    err(i) = norm( [chs(1) chs(2) 1]' - [UV(1,i) UV(2,i) 1]' );
end
RMS = RMSE_calc(err);