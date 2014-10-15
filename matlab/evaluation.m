function [err] = evaluation(H, H0, S, norm_xy)
% HCS-H0CS
n = size(S,3);
err = zeros(1,n);
%H_inv = inv(H);
%H0_inv = inv(H0);
for i=1:n
    s = S{i};
%     hs = H_inv' * (s/H);
%     h0s = H0_inv' * (s/H0);
%     chs = getEllipseCenter(hs);
%     ch0s = getEllipseCenter(h0s);
%     err(i) = norm(chs-ch0s);
    cs = getEllipseCenter(s);
    hcs = reshape(H, 3, 3)*[cs(1) cs(2) 1]';
    hcs = hcs / hcs(3);
    h0cs = H0*[cs(1) cs(2) 1]';
    h0cs = h0cs / h0cs(3);
    err(i) = norm( [hcs(1)/norm_xy hcs(2)/norm_xy 1]' - ...
       [h0cs(1)/norm_xy h0cs(2)/norm_xy 1]' );
end
