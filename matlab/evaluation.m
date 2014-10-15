% Copyright (c) <2014>, <vicrucann@gmail.com>
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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
