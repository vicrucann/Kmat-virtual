function [err] = geometrical_error(H, S, UV, norm_xy)
% HCS-CH0S(=UV)
n = size(S,3);
err = zeros(1,n);
for i = 1 : n
    s = S{i};
    cs = getEllipseCenter(s);
    hcs = reshape(H, 3, 3)*[cs(1) cs(2) 1]'; %[CS(1,i)/norm_xy CS(2,i)/norm_xy 1]';
    hcs = hcs/hcs(3);
    err(i) = sqrt((hcs(1)/norm_xy-UV(1,i)/norm_xy)^2 + (hcs(2)/norm_xy-UV(2,i)/norm_xy)^2);
    %err(i) = norm([hcs(1)/norm_xy hcs(2)/norm_xy 1]' - [UV(1,i)/norm_xy UV(2,i)/norm_xy 1]');
end