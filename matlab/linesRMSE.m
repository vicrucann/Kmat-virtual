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

function [RMS RMS_max flag_in] = linesRMSE(points, dx, dy, w, h, plotfigure)
if dy>dx
    UV = zeros(2,dx*dy);
    for i = dx:-1:1
        line = points(:,i:dx:dx*dy-dx+i);
        UV(:,(dx-i)*dy+1:(dx-i)*dy+dy) = line;
    end
    points = UV;
    tmp = dx;
    dx = dy;
    dy = tmp;
end

di = abs(dx-dy)+1;
% all the horizontal lines
rms_hor = zeros(dx,dy);
rmse_hor = 0;
step = 25;
if plotfigure == 1
    figure;
    hold on;
end
for y = 1:dy
    line = points(:,(y-1)*dx+1:dx*y);
    rms_hor(:,y) = residuals(line,w,h);
    if plotfigure == 1
        plot(line(1,:), line(2,:), 'Color', 'blue');
        for k = 1:dx
            if abs(rms_hor(k,y)) > 0.1
                text(line(1,k), line(2,k)+step, num2str(rms_hor(k,y)), 'Color', 'blue', 'FontSize', 7);
            end
        end
    end
    rmse_hor = rmse_hor + rms_hor(:,y)' * rms_hor(:,y);
end
rmse_hor_max = max(max(abs(rms_hor)));
rmse_hor = sqrt(rmse_hor/(dy*dx));

% all the vertical lines
rms_ver = zeros(dy,dx);
rmse_ver = 0;
for x = 1:dx
    line = points(:,x:dx:dx*dy-dx+x);
    rms_ver(:,x) = residuals(line,w,h);
    if plotfigure == 1
        plot(line(1,:), line(2,:), 'Color', 'red');
        for k = 1:dy
            if abs(rms_ver(k,x)) > 0.1
                text(line(1,k), line(2,k)-step, num2str(rms_ver(k,x)), 'Color', 'red', 'FontSize', 7);
            end
        end
    end
    rmse_ver = rmse_ver + rms_ver(:,x)' * rms_ver(:,x);
end
rmse_ver_max = max(max(abs(rms_ver)));
rmse_ver = sqrt(rmse_ver/(dy*dx));

% all diagonal lines from left to right
rms_di1 = zeros(dy,di);
rmse_di1 = 0;
for x = 1:di
    line = points(:,x:dx+1:dx*dy-di+x);
    rms_di1(:,x) = residuals(line,w,h);
    if plotfigure == 1
        plot(line(1,:), line(2,:), 'Color', 'green');
        for k = 1:dy
            if abs(rms_di1(k,x)) > 0.1
                text(line(1,k)-step*4, line(2,k), num2str(rms_di1(k,x)), 'Color', 'green', 'FontSize', 7);
            end
        end
    end
    rmse_di1 = rmse_di1 + rms_di1(:,x)' * rms_di1(:,x);
end
rmse_di1_max = max(max(abs(rms_di1)));
rmse_di1 = sqrt(rmse_di1/(dy*di));

% all diagonal lines from right to left
rms_di2 = zeros(dy,di);
rmse_di2 = 0;
for x = dx:-1:dx-di+1
    line = points(:,x:dx-1:dx*dy-dx+di);
    rms_di2(:,x) = residuals(line,w,h);
    if plotfigure == 1
        plot(line(1,:), line(2,:), 'Color', 'magenta');
        for k = 1:dy
            if abs(rms_di2(k,x))  > 0.1
                text(line(1,k), line(2,k), num2str(rms_di2(k,x)), 'Color', 'magenta', 'FontSize', 7);
            end
        end
    end
    rmse_di2 = rmse_di2 + rms_di2(:,x)' * rms_di2(:,x);
end
rmse_di2_max = max(max(abs(rms_di2)));
rmse_di2 = sqrt(rmse_di2/(dy*di));

if plotfigure == 1
    title('distortion correction error');
    hold off;
end

RMS = mean([rmse_hor rmse_ver rmse_di1 rmse_di2]);
RMS_max = max([rmse_hor_max rmse_ver_max rmse_di1_max rmse_di2_max]);

% now find and eliminate ourliers
thre_out = RMS_max; %(RMS_max-RMS)/4+RMS;
flag_in = ones(1,dx*dy); % 0 for ourlier, otherwise 1
% among horizontal lines
for y = 1:dy
    line = points(:,(y-1)*dx+1:dx*y);
    res = residuals(line,w,h);
    idx0 = find(abs(res)>thre_out);
    for k = 1:size(idx0,1)
        flag_in(1,idx0(k)+(y-1)*dx) = 0;
    end
end

% among vertical lines
for x = 1:dx
    line = points(:,x:dx:dx*dy-dx+x);
    res = residuals(line,w,h);
    idx0 = find(abs(res)>thre_out);
    for k = 1:size(idx0,1)
        flag_in(1,x+(idx0(k)-1)*dx) = 0;
    end
end

% among diagonals from left to right
for x = 1:di
    line = points(:,x:dx+1:dx*dy-di+x);
    res = residuals(line,w,h);
    idx0 = find(abs(res)>thre_out);
    for k = 1:size(idx0,1)
        flag_in(1,dx-idx0(k)+1) = 0;
    end
end

% among diagonals from left to right
for x = dx:-1:dx-di+1
    line = points(:,x:dx-1:dx*dy-dx+di);
    res = residuals(line,w,h);
    idx0 = find(abs(res)>thre_out);
    for k = 1:size(idx0,1)
        flag_in(1,(idx0(k)-1)*(dx-1)+x) = 0;
    end
end

end

function RMS = residuals(line, w, h)
order = 3;
nPoints = size(line,2);
sizexy = (order+1)*(order+2)/2;
paramsX = zeros(1,sizexy); paramsX(1,sizexy-2) = 1;
paramsY = zeros(1,sizexy); paramsY(1,sizexy-1) = 1;
xp = double(w/2)+0.2;
yp = double(h/2)+0.2;
[coefTermX coefTermY] = coefTermCalc(line(1,:), line(2,:), order, xp, yp);
ux = bicubicDistModel(paramsX, coefTermX);
uy = bicubicDistModel(paramsY, coefTermY);

Ax = mean(ux);
Ay = mean(uy);
x = ux-Ax;
y = uy-Ay;
Axy = ux'*uy/nPoints;
Vxx = (ux-Ax)'*(ux-Ax)/nPoints; %norm(ux-Ax);
Vyy = (uy-Ay)'*(uy-Ay)/nPoints; %norm(uy-Ay);
Vxy = Axy-Ax*Ay;
Vxx_yy = Vxx-Vyy;
theta = 0.5*atan2(-2*Vxy, Vxx-Vyy);

RMS = sin(theta)*x+cos(theta)*y; %sqrt(0.5*(Vxx + Vyy - sqrt(Vxx_yy^2 + 4*Vxy^2) ));
end

function [coefTermX coefTermY] = coefTermCalc(pointX, pointY, deg, xp, yp)

nPoints = size(pointX,2);
scale = 0;
sizex = (deg + 1) * (deg + 2) / 2;
sizey = (deg + 1) * (deg + 2) / 2;
coefTermX = zeros(sizex, nPoints);
coefTermY = zeros(sizey, nPoints);
idx = 1;
for ii = deg : -1: 0
	for j = 0 : ii
		for k = 1 : nPoints
			coefTermX(idx, k) = (pointX(1,k)-xp)^(ii-j) * (pointY(1,k)-yp)^j;
            coefTermY(idx, k) = coefTermX(idx, k);
		end
		idx = idx + 1;
	end
end

end

function model = bicubicDistModel(param, coefTerm)
model = coefTerm'*param';
end
