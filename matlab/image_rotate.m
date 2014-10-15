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

function [I H0S CH0S H0CS] = image_rotate(I0, S, H0, pix_flag, wim, hem)
I = 0;
if pix_flag == 1
    [he wi] = size(I0);
    %wim = wi; hem = he;
    [X,Y]=meshgrid(1:wi,1:he);
    [Xi,Yi]=meshgrid(1:wim,1:hem);
    for x = 1:wim
        for y = 1:hem
            m = [x y 1]';
            M = H0\m;
            M = M/M(3);
            if (M(1)<1); M(1)=1; end;
            if (M(2)<1); M(2)=1; end;
            if (M(1)>wi); M(1)=wi; end;
            if (M(2)>he); M(2)=he; end;
            Xi(y,x) = M(1);
            Yi(y,x) = M(2);
        end
    end
    I = uint8(interp2(X,Y,I0,Xi,Yi,'cubic'));
end

n = size(S,3);
H0S = cell(3,3,n);
CH0S = zeros(2,n);
H0CS = zeros(2,n);
for i = 1:n
    s = S{i};
    cs = getEllipseCenter(s);
    h0s = (inv(H0))' * (s / H0);
    H0S{i} = h0s;
    CH0S(:,i) = getEllipseCenter(h0s);
    xy = H0 * [cs(1) cs(2) 1]';
    H0CS(:,i) = [xy(1)/xy(3) xy(2)/xy(3)];
end
