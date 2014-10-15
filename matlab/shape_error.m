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

function [err jac] = shape_error(H, S, UV, norm_xy)
% CHS-CH0S(=UV)
n = size(S,3);
err = zeros(1,n);
for i = 1 : n
    s = S{i};
    hs = H'*s*H;
    chs = getEllipseCenter(hs) / norm_xy;
    err(i) = sqrt((chs(1,1)-UV(1,i)/norm_xy)^2 + (chs(2,1)-UV(2,i)/norm_xy)^2);
end
jac = zeros(n,9);
for i = 1:n
    jac(i,:) = 0;%derivatives(H, S{i}, UV(1,i), UV(2,i));
end

function dEdh = derivatives(H, S, u, v)
h1 = H(1,1); h2 = H(1,2); h3 = H(1,3);
h4 = H(2,1); h5 = H(2,2); h6 = H(2,3);
h7 = H(3,1); h8 = H(3,2); h9 = H(3,3);
a = S(1,3); b = S(2,3); f = S(3,3);
dEdh = zeros(1,9);

dEdh(1,1) = (2*(-(((h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))*(h3 + a*h9))/...
    (-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + ...
    (h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) + ...
    ((h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))*(-((h2 + a*h8)*(h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)) - ...
    (h2 + a*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8)) + (2*h1 + 2*a*h7)*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))*...
    (h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/ Power(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*...
    (h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + (h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*...
    (h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)),2))*(-(((h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))*...
    (h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/...
    (-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*...
    (h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + ...
    (h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*...
    (h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) - u) + ...
    2*(-(((-(h1*(h2 + a*h8)) - h4*(h5 + b*h8) - h7*(a*h2 + b*h5 + f*h8))*(h3 + a*h9))/...
    (-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*...
    (h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + ...
    (h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*...
    (h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) + ...
    ((-(h1*(h2 + a*h8)) - h4*(h5 + b*h8) - h7*(a*h2 + b*h5 + f*h8))*...
    (-((h2 + a*h8)*(h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)) - ...
    (h2 + a*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8)) + ...
    (2*h1 + 2*a*h7)*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))*...
    (h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/...
    Power(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*...
    (h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + ...
    (h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*...
    (h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)),2) - ...
    ((-h2 - a*h8)*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/...
    (-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*...
    (h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + ...
    (h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*...
    (h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))))* (-(((-(h1*(h2 + a*h8)) - h4*(h5 + b*h8) - h7*(a*h2 + b*h5 + f*h8))*...
    (h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/  (-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*...
    (h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + (h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*...
    (h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) - v))/(2.*Sqrt(Power(-(((h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))*...
    (h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*...
    (h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) +(h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*...
    (h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) - u,2) +Power(-(((-(h1*(h2 + a*h8)) - h4*(h5 + b*h8) - h7*(a*h2 + b*h5 + f*h8))*...
    (h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*...
    (h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + (h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*...
    (h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) - v,2)));

dEdh(1,2) = (2*(((h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))*...
    ((h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(2*h2 + 2*a*h8) -...
    (h1 + a*h7)*(h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8) - (h1 + a*h7)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))...
    )*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/...
    Power(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + ...
    (h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)),2) - ...
    ((2*h2 + 2*a*h8)*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/...
    (-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + ...
    (h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))))*...
    (-(((h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/...
    (-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + ...
    (h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) - u) + ...
    2*(((-(h1*(h2 + a*h8)) - h4*(h5 + b*h8) - h7*(a*h2 + b*h5 + f*h8))*...
    ((h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(2*h2 + 2*a*h8) - ...
    (h1 + a*h7)*(h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8) - (h1 + a*h7)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))...
    )*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/...
    Power(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + ...
    (h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)),2) -...
    ((-h1 - a*h7)*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/...
    (-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + ...
    (h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))))*...
    (-(((-(h1*(h2 + a*h8)) - h4*(h5 + b*h8) - h7*(a*h2 + b*h5 + f*h8))*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/...
    (-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + ...
    (h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) - v))/...
    (2.*Sqrt(Power(-(((h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/...
    (-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + ...
    (h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) - u,2) +...
    Power(-(((-(h1*(h2 + a*h8)) - h4*(h5 + b*h8) - h7*(a*h2 + b*h5 + f*h8))*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/...
    (-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + ...
    (h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) - v,2)));

dEdh(1,3) = ((-2*(h1 + a*h7)*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))*...
    (-(((h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/...
    (-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + ...
    (h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) - u))/...
    (-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + ...
    (h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))) - ...
    (2*(h1 + a*h7)*(-(h1*(h2 + a*h8)) - h4*(h5 + b*h8) - h7*(a*h2 + b*h5 + f*h8))*...
    (-(((-(h1*(h2 + a*h8)) - h4*(h5 + b*h8) - h7*(a*h2 + b*h5 + f*h8))*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/...
    (-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + ...
    (h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) - v))/...
    (-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + ...
    (h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))))/...
    (2.*Sqrt(Power(-(((h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/...
    (-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + ...
    (h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) - u,2) + ...
    Power(-(((-(h1*(h2 + a*h8)) - h4*(h5 + b*h8) - h7*(a*h2 + b*h5 + f*h8))*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/...
    (-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + ...
    (h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) - v,2)));

dEdh(1,4) = (2*(-(((h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))*(h6 + b*h9))/...
    (-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + ...
    (h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) + ...
    ((h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))*...
    (-((h5 + b*h8)*(h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)) - ...
    (h5 + b*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8)) + ...
    (2*h4 + 2*b*h7)*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/...
    Power(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + ...
    (h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)),2))*...
    (-(((h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/...
    (-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + ...
    (h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) - u) +...
    2*(-(((-(h1*(h2 + a*h8)) - h4*(h5 + b*h8) - h7*(a*h2 + b*h5 + f*h8))*(h6 + b*h9))/...
    (-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + ...
    (h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) + ...
    ((-(h1*(h2 + a*h8)) - h4*(h5 + b*h8) - h7*(a*h2 + b*h5 + f*h8))*...
    (-((h5 + b*h8)*(h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)) - ...
    (h5 + b*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8)) + ...
    (2*h4 + 2*b*h7)*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/...
    Power(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + ...
    (h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)),2) - ...
    ((-h5 - b*h8)*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/...
    (-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + ...
    (h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))))*...
    (-(((-(h1*(h2 + a*h8)) - h4*(h5 + b*h8) - h7*(a*h2 + b*h5 + f*h8))*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/...
    (-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + ...
    (h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) - v))/...
    (2.*Sqrt(Power(-(((h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/...
    (-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + ...
    (h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) - u,2) +...
    Power(-(((-(h1*(h2 + a*h8)) - h4*(h5 + b*h8) - h7*(a*h2 + b*h5 + f*h8))*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/...
    (-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + ...
    (h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) - v,2)));

dEdh(1,5) = (2*(((h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))*...
    ((h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(2*h5 + 2*b*h8) -...
    (h4 + b*h7)*(h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8) - (h4 + b*h7)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))...
    )*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/...
    Power(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + ...
    (h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)),2) - ...
    ((2*h5 + 2*b*h8)*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/...
    (-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + ...
    (h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))))*...
    (-(((h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/...
    (-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + ...
    (h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) - u) + ...
    2*(((-(h1*(h2 + a*h8)) - h4*(h5 + b*h8) - h7*(a*h2 + b*h5 + f*h8))*...
    ((h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(2*h5 + 2*b*h8) - ...
    (h4 + b*h7)*(h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8) - (h4 + b*h7)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))...
    )*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/...
    Power(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + ...
    (h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)),2) - ...
    ((-h4 - b*h7)*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/...
    (-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + ...
    (h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))))*...
    (-(((-(h1*(h2 + a*h8)) - h4*(h5 + b*h8) - h7*(a*h2 + b*h5 + f*h8))*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/...
    (-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + ...
    (h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) - v))/...
    (2.*Sqrt(Power(-(((h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/...
    (-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + ...
    (h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) - u,2) + ...
    Power(-(((-(h1*(h2 + a*h8)) - h4*(h5 + b*h8) - h7*(a*h2 + b*h5 + f*h8))*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/...
    (-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + ...
    (h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) - v,2)));

dEdh(1,6) = ((-2*(h4 + b*h7)*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))*...
    (-(((h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/...
    (-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + ...
    (h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) - u))/...
    (-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + ...
    (h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))) - ...
    (2*(h4 + b*h7)*(-(h1*(h2 + a*h8)) - h4*(h5 + b*h8) - h7*(a*h2 + b*h5 + f*h8))*...
    (-(((-(h1*(h2 + a*h8)) - h4*(h5 + b*h8) - h7*(a*h2 + b*h5 + f*h8))*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/...
    (-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + ...
    (h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) - v))/...
    (-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + ...
    (h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))))/...
    (2.*Sqrt(Power(-(((h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/...
    (-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + ...
    (h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) - u,2) + ...
    Power(-(((-(h1*(h2 + a*h8)) - h4*(h5 + b*h8) - h7*(a*h2 + b*h5 + f*h8))*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/...
    (-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + ...
    (h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) - v,2)));

dEdh(1,7) = (2*(-(((h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))*(a*h3 + b*h6 + f*h9))/...
    (-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + ...
    (h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) + ...
    ((h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))*...
    (-((a*h2 + b*h5 + f*h8)*(h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)) - ...
    (a*h2 + b*h5 + f*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8)) + ...
    (2*a*h1 + 2*b*h4 + 2*f*h7)*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))*...
    (h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/...
    Power(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + ...
    (h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)),2))*...
    (-(((h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/...
    (-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + ...
    (h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) - u) + ...
    2*(-(((-(h1*(h2 + a*h8)) - h4*(h5 + b*h8) - h7*(a*h2 + b*h5 + f*h8))*(a*h3 + b*h6 + f*h9))/...
    (-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + ...
    (h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) + ...
    ((-(h1*(h2 + a*h8)) - h4*(h5 + b*h8) - h7*(a*h2 + b*h5 + f*h8))*...
    (-((a*h2 + b*h5 + f*h8)*(h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)) - ...
    (a*h2 + b*h5 + f*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8)) + ...
    (2*a*h1 + 2*b*h4 + 2*f*h7)*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))*...
    (h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/...
    Power(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + ...
    (h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)),2) - ...
    ((-(a*h2) - b*h5 - f*h8)*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/...
    (-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + ...
    (h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))))*...
    (-(((-(h1*(h2 + a*h8)) - h4*(h5 + b*h8) - h7*(a*h2 + b*h5 + f*h8))*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/...
    (-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + ...
    (h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) - v))/...
    (2.*Sqrt(Power(-(((h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/...
    (-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + ...
    (h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) - u,2) + ...
    Power(-(((-(h1*(h2 + a*h8)) - h4*(h5 + b*h8) - h7*(a*h2 + b*h5 + f*h8))*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/...
    (-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + ...
    (h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) - v,2)));

dEdh(1,8) = (2*(((h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))*...
    ((h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(2*a*h2 + 2*b*h5 + 2*f*h8) - ...
    (a*h1 + b*h4 + f*h7)*(h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8) - ...
    (a*h1 + b*h4 + f*h7)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8)))*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/...
    Power(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + ...
    (h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)),2) - ...
    ((2*a*h2 + 2*b*h5 + 2*f*h8)*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/...
    (-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + ...
    (h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))))*...
    (-(((h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/...
    (-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + ...
    (h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) - u) + ...
    2*(((-(h1*(h2 + a*h8)) - h4*(h5 + b*h8) - h7*(a*h2 + b*h5 + f*h8))*...
    ((h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(2*a*h2 + 2*b*h5 + 2*f*h8) - ...
    (a*h1 + b*h4 + f*h7)*(h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8) - ...
    (a*h1 + b*h4 + f*h7)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8)))*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/...
    Power(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + ...
    (h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)),2) - ...
    ((-(a*h1) - b*h4 - f*h7)*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/...
    (-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + ...
    (h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))))*...
    (-(((-(h1*(h2 + a*h8)) - h4*(h5 + b*h8) - h7*(a*h2 + b*h5 + f*h8))*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/...
    (-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + ...
    (h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) - v))/...
    (2.*Sqrt(Power(-(((h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/...
    (-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + ...
    (h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) - u,2) + ...
    Power(-(((-(h1*(h2 + a*h8)) - h4*(h5 + b*h8) - h7*(a*h2 + b*h5 + f*h8))*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/...
    (-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + ...
    (h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) - v,2)));

dEdh(1,9) = ((-2*(a*h1 + b*h4 + f*h7)*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))*...
    (-(((h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/...
    (-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + ...
    (h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) - u))/...
    (-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + ...
    (h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))) - ...
    (2*(a*h1 + b*h4 + f*h7)*(-(h1*(h2 + a*h8)) - h4*(h5 + b*h8) - h7*(a*h2 + b*h5 + f*h8))*...
    (-(((-(h1*(h2 + a*h8)) - h4*(h5 + b*h8) - h7*(a*h2 + b*h5 + f*h8))*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/...
    (-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + ...
    (h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) - v))/...
    (-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + ...
    (h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))))/...
    (2.*Sqrt(Power(-(((h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/...
    (-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + ...
    (h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) - u,2) + ...
    Power(-(((-(h1*(h2 + a*h8)) - h4*(h5 + b*h8) - h7*(a*h2 + b*h5 + f*h8))*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/...
    (-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + ...
    (h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) - v,2)));
          
function res = Sqrt(exp)
res = sqrt(exp);

function res = Power(base, pow)
res = base^pow;
