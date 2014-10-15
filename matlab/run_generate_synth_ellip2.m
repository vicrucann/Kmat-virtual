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

clc;
clear all;
close all;

% generate synthetic data using some synthetic ground truth
% do not include distortion, but rather noise level for coordinate
% use ellip method to obtain homography, compare with ground truth
% compare UV - CH0S - H0CS

%% PATTERN
scale = 2; % 1.5
trials = 10;
wi_pix = 3888/3; % real image size
he_pix = 2592/3;
WI_pix = wi_pix * scale; % for the pattern
HE_pix = he_pix * scale;

ncircleX = 14; ncircleY = 10;
ncircle = ncircleX*ncircleY;
dX = 3.0; dY = 3.0;
CM1 = 100;
MM1 = CM1/10;
WI = dX*ncircleX; HE = dY*ncircleY; % in cm
RADI = min(dX,dY)/3; % cm
begin = 0;
DELTA = 0; delta = DELTA/scale;
XC = 0.5*(dX^2-begin)/(dX-begin) + DELTA;
YC = 0.5*(dY^2-begin)/(dY-begin) + DELTA;
x_kpts = (0 : dX : WI-dX);
y_kpts = (0 : dY : HE-dY);

%scale_cm2pix = floor(min(WI_pix/(CM1*WI), HE_pix/(CM1*HE)))*CM1;
%scale_cm2pix = 150;
scale_cm2pix = 55;
begin_pix_x = 0.5*(WI_pix - WI*scale_cm2pix);
begin_pix_y = 0.5*(HE_pix - HE*scale_cm2pix);
dX_pix = dX*scale_cm2pix;
dY_pix = dY*scale_cm2pix;
begin_pix = begin*scale_cm2pix; % 1:w, 1:h
xc = XC * scale_cm2pix;
yc = YC * scale_cm2pix;
DELTA_pix = DELTA*scale_cm2pix;
x_kpts_pix = x_kpts * scale_cm2pix + begin_pix_x;
y_kpts_pix = y_kpts * scale_cm2pix + begin_pix_y;

if exist('I0.pgm') == 0
    fprintf('Building an image...\n');
    I0 = zeros(HE_pix, WI_pix);
    I0(:,:) = 255;
end
ground = zeros(2,ncircle);
idx = 1;
S = cell(3,3,ncircle);
Sn = cell(3,3,ncircle);
Sc = cell(3,3,4);
N1 = [2/WI 0 -1; 0 2/HE -1; 0 0 1];
idx_corner = 1;
idx_i = 1;
idx_j = 1;
for j = y_kpts_pix
    idx_i = 1;
    for i = x_kpts_pix
        XCi = XC + x_kpts(idx_i)-begin; 
        YCi = YC + y_kpts(idx_j)-begin;
        xci = xc + i-begin_pix;
        yci = yc + j-begin_pix;
        ground(:,idx) = [xci yci]';
        pn = N1*[XCi; YCi; 1]; % cm
        S{idx} = getCircleMatrix(XCi, YCi, RADI); % cm
        Sn{idx} = getCircleMatrix(pn(1), pn(2), 2*RADI/min(WI,HE)); % cm
        if exist('I0.pgm') == 0
            for y = j : j+dY_pix-1 % pix
                for x = i : i+dX_pix-1
                    if ((x-xci)^2 + (y-yci)^2 <= (RADI*scale_cm2pix)^2)
                        I0(y, x) = 0;
                    end
                end
            end
        end
        if (i == x_kpts_pix(1) && j == y_kpts_pix(1)) || (i == x_kpts_pix(end) && j == y_kpts_pix(1)) ... % collect all the corner data
                || (i == x_kpts_pix(1) && j == y_kpts_pix(end)) || (i == x_kpts_pix(end) && j == y_kpts_pix(end))
            Sc{idx_corner} = S{idx};
            idx_corner = idx_corner + 1;
        end
        idx = idx+1;
        idx_i = idx_i+1;
    end
    idx_j = idx_j+1;
end
% save pattern image
if exist('I0.pgm') == 0
    imwrite(uint8(I0), 'I0.pgm', 'pgm');
else
    I0 = imread('I0.pgm');
end

%% SYNTHETIC DISTORTION-FREE DATA
nset = 5;
nimage = 5;

[ny nx] = size(I0);
[mx my] = meshgrid(1:nx, 1:ny);
[nnx,nny] = size(mx);
px = reshape(mx,nnx*nny,1);
py = reshape(my,nnx*nny,1);

scale_k = 1;
alpha = 1250*scale_k; beta = 1250*scale_k; gamma = 1.09083; % K matrix params
u0 = 0.5*wi_pix; v0 = 0.5*he_pix;

K0 = [alpha gamma u0; ...
    0 beta v0; ...
    0 0 1];

Thetas = zeros(3, nimage); % camera rotations
Thetas(:, 1) = [deg2rad(3) deg2rad(2) deg2rad(1)]';
Thetas(:, 2) = -[deg2rad(0.5) deg2rad(25) deg2rad(0.1)]';
Thetas(:, 3) = [deg2rad(0.7) deg2rad(20) deg2rad(0.4)]';%1/sqrt(5)*[deg2rad(-30) deg2rad(-30) deg2rad(-15)]';
Thetas(:, 4) = [deg2rad(30) deg2rad(2) deg2rad(0.6)]';
Thetas(:, 5) = -[deg2rad(25) deg2rad(7) deg2rad(0.5)]';
Ts = zeros(3, nimage); % camera translations
Ts(:, 1) = scale_k*[-WI/2      -HE/2       WI*alpha/wi_pix+8]'; % in CM
Ts(:, 2) = scale_k*[-WI/2      -HE/2       WI*alpha/wi_pix+10]';
Ts(:, 3) = scale_k*[-WI/2-2    -HE/2       WI*alpha/wi_pix+20]';
Ts(:, 4) = scale_k*[-WI/2-2    -HE/2       WI*alpha/wi_pix+8]';
Ts(:, 5) = scale_k*[-WI/2      -HE/2       WI*alpha/wi_pix+15]';
theta_noise = [deg2rad(4) deg2rad(2) deg2rad(0.8)]';
trans_noise = [0.5 0.5 1]'; % in cm


H0 = cell(3,3,nset*nimage);
CH0S = cell(2,ncircle,nset*nimage);
H0CS = cell(2,ncircle,nset*nimage);
kernel_avg = ones(scale, scale)/scale^2;
for s = 1:nset
    %var = load(['images\Hellip_set' num2str(s) '.mat']);
    %Hs = var.H_ellip;
    for i = 1:nimage
        rand_rat_th = zeros(3,1);
        rand_rat_tr = zeros(3,1);
        for n=1:3
            rand_rat_th(n) = rand(1,1)*2-1;
            rand_rat_tr(n) = rand(1,1)*2-1;
        end
        theta = Thetas(:,i) + rand_rat_th.*theta_noise;
        trans = Ts(:,i) + rand_rat_tr.*trans_noise;
        h0 = groundH(alpha, beta, gamma, u0, v0, theta, trans);
        if exist(['images\H' num2str(s) num2str(i) '.mat']) == 0
            save(['images\H' num2str(s) num2str(i) '.mat'], 'h0');
        end
        var = load(['images\H' num2str(s) num2str(i) '.mat']);
        h0 = var.h0;
        H0{(s-1)*nimage+i} = h0;
        ch0s = zeros(2,ncircle);
        h0cs = zeros(2,ncircle);
        for t = 1:ncircle
            ch0s(:,t) = getEllipseCenter( inv(h0)' * S{t} / h0 );
            cs = getEllipseCenter(S{t});
            xy = h0 * [cs(1) cs(2) 1]';
            h0cs(:,t) = [xy(1)/xy(3) xy(2)/xy(3)];
        end
        CH0S{(s-1)*nimage+i} = ch0s;
        H0CS{(s-1)*nimage+i} = h0cs;
        if exist(['images\Iu' num2str(s) num2str(i) '.pgm']) == 0
            idx_interp = 1;
            % in order to interpolate the right values, we go backwards
            % big_img-> small_img-> undistort,inv(H)-> pattern,cm-> big_img
            fprintf('set %i image %i...', s, i);
            p_pix_big = [px py]';
            p_pix_small = p_pix_big./scale;
   
            p_cm = h0 \ [p_pix_small; ones(1,length(p_pix_small))];
            p_cm = [p_cm(1,:)./p_cm(3,:);  p_cm(2,:)./p_cm(3,:)];
            p_pix_big_col = p_cm(1:2,:) * scale_cm2pix;
            px_interp = p_pix_big_col(1,:) + begin_pix_x;
            py_interp = p_pix_big_col(2,:) + begin_pix_y;
            
            % interpolate 
            mx_interp = reshape(px_interp, [ny nx]); 
            my_interp = reshape(py_interp, [ny nx]);
            Iu = interp2(mx, my, double(I0), mx_interp, my_interp);
            imwrite(uint8(Iu), ['images\Iu' num2str(s) num2str(i) '.pgm'], 'pgm');
            fprintf('done.\n'); 
        end
    end
    
    % smooth and scale down
    for i = 1:nimage
    if exist(['images\_iu' num2str(s) num2str(i) '.pgm']) == 0
        Iu = imread(['images\Iu' num2str(s) num2str(i) '.pgm']);
        Iu_avg = conv2(double(Iu), kernel_avg, 'same');
        iu = imresize(uint8(Iu_avg), 1/scale);
        %iu = uint8(Iu_avg(scale/2 : scale : end, scale/2 : scale : end));
        imwrite(iu, ['images\_iu' num2str(s) num2str(i) '.pgm'], 'pgm');
    end
    end
end

%% DETECT KEYPOINTS 'UV', COMAPARE THEM WITH GROUNDTRUTH 'CH0S', +'H0CS'
for s = 1:nset
    for i = 1:nimage
        if exist(['images\_UV' num2str(s) num2str(i) '.mat'])==0
            Ii = imread(['images\_iu' num2str(s) num2str(i) '.pgm']);
            UV = detect_centers(Ii, S, Sc);
            save(['images\_UV' num2str(s) num2str(i) '.mat'], 'UV');
        end
    end
end
fprintf('\n KEYPOINT DETECTION RESULTS \n');
for s = 1:nset
    for i = 1:nimage
        ch0s = CH0S{(s-1)*nimage+i};
        h0cs = H0CS{(s-1)*nimage+i};
        var = load(['images\_UV' num2str(s) num2str(i) '.mat']);
        UV = var.UV;
        fprintf('Set %i Img %i detected points max(CH0S-UV):  dx = %f,    dy = %f\n', s, i,...
            max(abs(ch0s(1,:)-UV(1,:))), max(abs(ch0s(2,:)-UV(2,:))));
        %fprintf('Set %i Img %i detected points max(H0CS-UV):  dx = %f,    dy = %f\n', s, i,...
        %    max(abs(h0cs(1,:)-UV(1,:))), max(abs(h0cs(2,:)-UV(2,:))));
    end
end

%% EXTRACT HOMOGRAPHIES
fprintf('\n HOMOGRAPHY RESULTS \n');
options = optimset('Algorithm', 'levenberg-marquardt', 'ScaleProblem', 'Jacobian', 'Display', 'off', ...
    'Jacobian', 'off', 'FinDiffType', 'central', 'DerivativeCheck', 'off');
lb = []; ub = [];
norm_xy = 1;
nnoise = 11; noise_dev_max = 0.1; noise_dev = (0: noise_dev_max/(nnoise-1) :noise_dev_max);
%H_ellip = cell(3,3,nimage*nset);
H_ellip = zeros(3,3,nnoise,nimage,nset);
for s = 1:nset
    for i = 1:nimage
        var = load(['images\_UV' num2str(s) num2str(i) '.mat']);
        UV = var.UV;
        %UV = CH0S{(s-1)*nimage+i};
        h0 = H0{(s-1)*nimage+i};
        for n = 1:nnoise
            noise_int = [-noise_dev(n) noise_dev(n)];
            noise_val = randn(size(UV));
            noise_val = noise_int(1)+diff(noise_int)*(noise_val-min(noise_val(:)))/(max(noise_val(:))-min(noise_val(:)));
            UV_noise = UV + noise_val;
            
            h_ini = homography_solve(S, UV_noise);
            h_ini_inv = inv(h_ini);
            h_fin_ellip_inv = lsqnonlin(@(x) shape_error(x,S, UV_noise, norm_xy), h_ini_inv, lb, ub, options);
            h_ellip = inv(h_fin_ellip_inv);
            H_ellip(:,:,n,i,s) = h_ellip;
            %H_ellip{(s-1)*nimage+i} = h_ellip;    
        end
        RMSE_ellip_ini = RMSE_calc(evaluation(h_ini, h0, S, norm_xy));
        RMSE_ellip_fin = RMSE_calc(evaluation(h_ellip, h0, S, norm_xy));
        RMS_ellip = calc_reprojection(UV, h_ellip, S);
        fprintf('Set %i img %i |CHS-CH0S| RMSE_{ini} = %f,   RMSE_{estim} = %f\n',s,i,...
            RMSE_ellip_ini, RMSE_ellip_fin);
        fprintf('Reprojection |UV-CHS| RMS: %f\n', RMS_ellip);
    end
end

%% EXTRACT K MATRIX, GET STATS
fprintf('\n K MATRIX RESULTS \n');
alpha_ellip = zeros(nset,1); beta_ellip = zeros(nset,1);
gamma_ellip = zeros(nset,1); u0_ellip = zeros(nset,1); v0_ellip = zeros(nset,1);
alpha_std = zeros(1,nnoise); beta_std = zeros(1,nnoise);
gamma_std = zeros(1,nnoise); u0_std = zeros(1,nnoise); v0_std = zeros(1,nnoise);
for n = 1:nnoise
    for s = 1:nset
        H_set = cell(3,3,nimage);
        for i = 1:nimage
            %H_set{i} = H0{(s-1)*nimage+i};
            H_set{i} = H_ellip(:,:,n,i,s);%H_ellip{(s-1)*nimage+i};
        end
        K_ellip = extractKfromH(H_set);
        [eae ebe eue eve ege] = calc_Kparams_errors(K0, K_ellip);
        fprintf('Set %i |K_{ellip}-K_{g.t}| errors:   %f   %f   %f   %f   %f\n', s, eae, ebe, eue, eve, ege);
        
        alpha_ellip(s) = K_ellip(1,1);
        beta_ellip(s) = K_ellip(2,2);
        gamma_ellip(s) = K_ellip(1,2);
        u0_ellip(s) = K_ellip(1,3);
        v0_ellip(s) = K_ellip(2,3);
    end
    alpha_std(1,n) = std(alpha_ellip);  beta_std(1,n) = std(beta_ellip);
    gamma_std(1,n) = std(rad2deg(acos(gamma_ellip./beta_ellip)));
    u0_std(1,n) = std(u0_ellip); v0_std(1,n) = std(v0_ellip);
end

figure;
subplot(2,1,1);
hold on;
plot(noise_dev, alpha_std, 'Color', 'red');
plot(noise_dev, beta_std, 'Color', 'magenta');
legend('dev_{alpha}', 'dev_{beta}');
xlabel('noise level (pixels)'); ylabel('deviation of parameter');
title('Deviation vs. noise for ellip method');
hold off;
subplot(2,1,2);
hold on;
plot(noise_dev, u0_std, 'Color', 'red');
plot(noise_dev, v0_std, 'Color', 'magenta');
legend('dev_{u0}', 'dev_{v0}');
xlabel('noise level (pixels)'); ylabel('deviation of parameter');
%title('Deviation vs. noise for ellip method');
hold off;



%fprintf('\n deviations for %i sets \ndev(alpha) =   %f \ndev(beta) =    %f \ndev(gamma) =  %f \ndev(u0) =  %f \ndev(v0) =  %f\n', ...
%    nset, alpha_std, beta_std, gamma_std, u0_std, v0_std);






















