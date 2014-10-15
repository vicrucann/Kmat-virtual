function C = detect_centers(I1, S, Sc)
[HE WI] = size(I1);
%max_col = max(max(I1(I1~=255))); min_col = min(min(I1(I1~=0)));
%thre = double((max_col-min_col))*0.6/255;% 0.55;
thre = 0.55;
I1_bi = im2bw(I1, thre);
I1_bi = imfill(I1_bi,[1 1]);
%I2 = imcomplement(I1);
%I_bi = im2bw(I2, thre);
%I_bi = imfill(I_bi, 'holes');
I_bi = imfill(imcomplement(I1_bi), 'holes');
I2 = imcomplement(I1_bi);

I_la = logical(I_bi);
taches = regionprops(I_la, I2, 'all');
ntache = size(taches, 1);
radis = zeros(1, ntache);
center = zeros(2, ntache);
compactness = zeros(1,ntache);
% initilization
for k = 1 : ntache
    c = taches(k).Centroid;
    a = taches(k).MajorAxisLength/2;
    radis(1, k) = a;
    center(:,k) = [c(1) c(2)]';
    compactness(1,k) = 4*pi*taches(k).FilledArea / (taches(k).Perimeter)^2;
end
% tache_noisy_idx = find(radis<30);
% taches(tache_noisy_idx) = [];
% center(:, tache_noisy_idx) = [];
% radis(tache_noisy_idx) = [];
% compactness(tache_noisy_idx) = [];
% 
tache_noisy_idx = find(compactness < 0.6 | compactness > 1.4);
taches(tache_noisy_idx) = [];
center(:, tache_noisy_idx) = [];
radis(tache_noisy_idx) = [];
ntache = size(taches, 1);

fprintf('detected %i shapes. \n', ntache);
corners = detect_corners(center);
h_ini = homography_solve(Sc, corners);
[center radis] = sort_keypnts(h_ini, S, center, radis);

fprintf('centroids initialization done. \n');
% calculate LMA centers
C = zeros(2,ntache);
for k = 1 : ntache
    % create a sub-image, save it to pgm file
    radik = radis(k);
    sizek = round(radik / 0.4); % because R = 0.8*W/2
    x1 = round(center(1,k)-sizek/2);
    x2 = round(center(1,k)+sizek/2);
    y1 = round(center(2,k)-sizek/2);
    y2 = round(center(2,k)+sizek/2);
    if (x1 < 1); x1 = 1; end;
    if (x2 < 1); x2 = 1; end;
    if (y1 < 1); y1 = 1; end;
    if (y2 < 1); y2 = 1; end;
    if (x1 > WI); x1 = WI; end;
    if (x2 > WI); x2 = WI; end;
    if (y1 > HE); y1 = HE; end;
    if (y2 > HE); y2 = HE; end;
    Ik = I1(y1:y2, x1:x2);
    %fprintf('Press!\n'); figure(1); imshow(Ik); waitforbuttonpress;
    imwrite(uint8(Ik), 'test2.pgm', 'pgm');
    %imwrite(uint8(Ik), ['figures\ellipse_3\test' num2str(setn) num2str(imn) num2str(k) '.pgm'], 'pgm');
    [~,res] = system('exe\ellipse_params test2.pgm ');
    Pr = strread(res);
    centerk = [Pr(4) Pr(5)]; 
    C(:,k) = [x1+centerk(1), y1+centerk(2)]';
    percent = 100*k/ntache;
    if mod(k, round(0.2*ntache)) == 0
        fprintf('%i', percent);
    elseif mod(k, round(0.04*ntache)) == 0
        fprintf('.');
    end
end

fprintf('\nLMA centers done.\n');

function [c r] = sort_keypnts(H, S, center, radis)
n = size(center,2);
c = zeros(size(center));
r = zeros(1,n);
for i = 1:n
    s = S{i};
    cs = getEllipseCenter(s);
    idx_sort = find_keypnt(H,cs,center);
    c(:,i) = center(:,idx_sort);
    r(1,i) = radis(1,idx_sort);
end

function idx_min = find_keypnt(H, xy, center)
idx_min = 0;
M = [xy(1) xy(2) 1]';
m = H*M;
m = m/m(3);
minlen = 1000;
for i = 1 : size(center,2)
    len = euc_dist_points(m(1), m(2), center(1,i), center(2,i));
    if (len < minlen)
        minlen = len;
        idx_min = i;
    end
end

function corners = detect_corners(center)
c = center;
n = size(center, 2);
corner_idx = zeros(1,4);
ncorners = 4;
idxc = 1;
for i = 1:n
    corner_flags = zeros(1,n);
    [neigh1 neigh2] = nn2(c, i);
    for j = 1 : n
        if (j ~= neigh1 && j~= neigh2 && j ~= i)
            corner_flags(1,j) = corner_check(c,i,neigh1,neigh2,j,deg2rad(5));
        end
    end
    if (sum(corner_flags) == 0)
        ncorners = ncorners - 1;
        corner_idx(idxc) = i;
        idxc = idxc + 1;
    end
end
if (ncorners ~= 0)
    fprintf('more than 4 corners were detected\n');
end
corners = c(:,corner_idx);
[~,corner_idx] = sort(corners(2,:));
c = corners(:,corner_idx);
wi = 2; he = 2;
for j = 1 : he
    row = c(1, (j-1)*wi+1 : (j-1)*wi+wi);
    [~,rowidx] = sort(row);
    for i = 1:wi
        corners(:,(j-1)*wi+i) = c(:, rowidx(i)+(j-1)*wi);
    end
end

function inside = corner_check(c, i, neigh1, neigh2, idx, tol)
inside = 1;
if (c(2,neigh1) > c(2,neigh2)) % neigh1 must be always smaller
    tmp = neigh2;
    neigh2 = neigh1; 
    neigh1 = tmp;
end
a = [c(1,i)-c(1,neigh1) c(2,i)-c(2,neigh1)]; % vectors
b = [c(1,i)-c(1,neigh2) c(2,i)-c(2,neigh2)];
h = [c(1,i)-c(1,idx) c(2,i)-c(2,idx)];
a = a/norm(a);
b = b/norm(b);
h = h/norm(h);
theta_ab = real(acos(dot(a,b))); %reference angle
theta_ab = wrapTo2Pi(theta_ab);
theta_ah = real(acos(dot(a,h)));
theta_ah = wrapTo2Pi(theta_ah);
theta_bh = real(acos(dot(b,h)));
theta_bh = wrapTo2Pi(theta_bh);
if (theta_ah < theta_ab+tol && theta_bh < theta_ab+tol && theta_ab < deg2rad(160)) 
    inside = 0;
end

function [neigh1 neigh2] = nn2(c, idx)
n = size(c,2);
sorted_neigh = zeros(1,n); % keep track on indexes, too
for i = 1:n
        sorted_neigh(1,i) = euc_dist(c, idx, i);
end
[~,idxn] = sort(sorted_neigh);
n1 = idxn(2);
n2 = idxn(3);
n3 = idxn(4);
x = [c(1,idx)-c(1,n1) c(2,idx)-c(2,n1)];%vectors
y = [c(1,idx)-c(1,n2) c(2,idx)-c(2,n2)];
z = [c(1,idx)-c(1,n3) c(2,idx)-c(2,n3)];
x = x/norm(x);
y = y/norm(y);
z = z/norm(z);
theta_n1n2 = acos(dot(x,y));
theta_n1n3 = acos(dot(x,z));
theta_n2n3 = acos(dot(y,z));
theta_n1n2 = wrapTo2Pi(theta_n1n2);
theta_n1n3 = wrapTo2Pi(theta_n1n3);
theta_n2n3 = wrapTo2Pi(theta_n2n3);
[~,idxn] = max([theta_n1n2, theta_n1n3, theta_n2n3]);
if (idxn == 1)
    neigh1 = n1;
    neigh2 = n2;
elseif (idxn == 2)
    neigh1 = n1;
    neigh2 = n3;
else
    neigh1 = n2;
    neigh2 = n3;
end

function len = euc_dist(c, idx1, idx2)
len = sqrt((c(1,idx1)-c(1,idx2))^2 + (c(2,idx1)-c(2,idx2))^2);

function len = euc_dist_points(x1, y1, x2, y2)
len = sqrt((x1-x2)^2 + (y1-y2)^2);
