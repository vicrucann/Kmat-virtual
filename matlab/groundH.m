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


function [H oc0] = groundH(alpha, beta, gamma, u0, v0, thetas, ts)
thetax = thetas(1);
thetay = thetas(2);
thetaz = thetas(3);

K = [...
    alpha gamma u0; ...
    0 beta v0; ...
    0 0 1];
Rx = [...
    1 0 0; ...
    0 cos(thetax) -sin(thetax); ...
    0 sin(thetax) cos(thetax)];
Ry = [... 
    cos(thetay) 0 sin(thetay); ...
    0 1 0; ...
    -sin(thetay) 0 cos(thetay)];
Rz = [...
    cos(thetaz) -sin(thetaz) 0; ...
    sin(thetaz) cos(thetaz) 0; ...
    0 0 1];
R = Rx*Ry*Rz;

oc0 = -R\ts; % optical center ground truth

H = K * [R(:,1) R(:,2) ts];
