%% Initialization

clear all
close all
load coastlines

set(0, 'defaultaxesfontsize', 19)
set(0, 'defaultlinelinewidth', 2)


%% Read Lyon data

load('../data/lyon_data.mat')
mmd_thr1 = 25;
mmd_thr2 = 1000;

d0   = datenum('01-Jan-2019');
dd0  = datenum('01-Jun-2019')-d0+1;
day0 = 206;

idx  = find(mmd > mmd_thr1 & mmd < mmd_thr2);
Nd   = length(idx);

d_f  = round(day(idx));

if(idx(1) == 1)
    d_i = round(day(idx(2:end)-1));
    d_i = cat(2,day0,d_i);
else
    d_i = round(day(idx-1));
end

d_x  = zeros(1,0);
for n = 1:Nd
d_x  = cat(2,d_x,d_i(n):d_f(n));
end
d_x  = unique(d_x)-dd0+1;
Nd   = length(d_x);


%% Files parameters

root    = '../data/trajectories_500m/';

fn      = dir([root,'*']);
Nf      = length(fn);

yy      = zeros(1,Nf-2);
mm      = zeros(1,Nf-2);
dd      = zeros(1,Nf-2);

Nh      = 121;
Nmm     = 12;


%% Histogram bins

dbin_l  = .5;

lat_bin =  -90:dbin_l:90;
Nlat    = length(lat_bin);
Nlat1   = Nlat-1;

lon_bin = -360:dbin_l:360;
Nlon    = length(lon_bin);
Nlon1   = Nlon-1;

dbin_z  = 100;
z_bin   = 0:dbin_z:30000;
Nz      = length(z_bin)-1;

dbin_p  = 10;
p_bin   = 0:dbin_p:1200;
Np      = length(p_bin)-1;


%% Memory allocation & reading

hist_xy = zeros(Nlat,Nlon,Nmm);
hist_x  = zeros(Nlon1,Nh,Nmm);
hist_y  = zeros(Nlat1,Nh,Nmm);
hist_z  = zeros(Nz,Nh,Nmm);
hist_p  = zeros(Np,Nh,Nmm);
norm    = zeros(1,Nmm);

for nf = 3:Nf
 


    
    vec  = fn(nf).name;
    n    = nf-2;
    
    yy(n) = str2num( vec(6:9) );
    mm(n) = str2num( vec(10:11));
    dd(n) = str2num( vec(12:13));
    
    d1    = num2str(mm(n));
    d2    = num2str(dd(n));
    d3    = num2str(yy(n));
    
    daten = datenum([d1,'-',d2,'-',d3])-d0+1;
    

    tmp   = readmatrix([root,fn(nf).name]);
    dat   = tmp(35:end,:);

    lat_h = dat(:,10);
    lon_h = dat(:,11);
    z_h   = dat(:,12);
    p_h   = dat(:,13);
    
    lon_h(lon_h > 0) = lon_h(lon_h > 0)-360;
    
    hist_xy(:,:,mm(n)) = hist_xy(:,:,mm(n)) + hist3([lat_h lon_h], {lat_bin, lon_bin});
        
    vec = 1-dat(:,9);
    for nh = 1:Nh
        
        idx = find(vec == nh);
        hist_z(:,nh,mm(n)) = hist_z(:,nh,mm(n)) + histcounts(z_h(idx),   z_bin')';
        hist_p(:,nh,mm(n)) = hist_p(:,nh,mm(n)) + histcounts(p_h(idx),   p_bin')';
        hist_x(:,nh,mm(n)) = hist_x(:,nh,mm(n)) + histcounts(lon_h(idx), lon_bin')';
        hist_y(:,nh,mm(n)) = hist_y(:,nh,mm(n)) + histcounts(lat_h(idx), lat_bin')';
        
    end

    norm(mm(n)) = norm(mm(n)) + 1;
    
    disp(nf)
end


%%

% v1  = [5:10 17:22];
% n1  = length(v1);
% v2  = [1:4 11:16 23:24];
% n2  = length(v2);

v1 = [1:4 11:12];
v2 = [5:10];
n1 = length(v1);
n2 = length(v2);

n1 = sum(norm(v1));
n2 = sum(norm(v2));

thr = 0;

mat1_xy = squeeze(sum(hist_xy(:,:,v1),3)); %Wet season
mat2_xy = squeeze(sum(hist_xy(:,:,v2),3)); %Dry season
mat3_xy = squeeze(sum(hist_xy(:,:,:),3)); % Total

mat1_x  = squeeze(sum(hist_x(:,:,v1),3));
mat2_x  = squeeze(sum(hist_x(:,:,v2),3));

mat1_y  = squeeze(sum(hist_y(:,:,v1),3));
mat2_y  = squeeze(sum(hist_y(:,:,v2),3));

mat1_p  = squeeze(sum(hist_p(:,:,v1),3));
mat2_p  = squeeze(sum(hist_p(:,:,v2),3));

mat1_z  = squeeze(sum(hist_z(:,:,v1),3));
mat2_z  = squeeze(sum(hist_z(:,:,v2),3));



%% Lat-lon total

close all
figure('Position', [1 1 800 400]);
mat3_xy(mat3_xy < thr) = NaN;
contourf(lon_bin, lat_bin, log( mat3_xy), 100, 'linestyle', 'none')
geoshow(coastlat,coastlon,'Color','k')
% ylim([15 45])
ylim([0 90])
% xlim([-180 -120])
xlim([-280 -100])
% caxis([3 7])
xlabel('Longitude (ºE)')
ylabel('Latitude (ºN)')
cb1 = colorbar;


%% Lat-lon by season


close all
figure('Position', [1 1 800 700]);
subplot(2,1,1)
mat1_xy(mat1_xy < thr) = NaN;
contourf(lon_bin, lat_bin, log( mat1_xy), 100, 'linestyle', 'none')
geoshow(coastlat,coastlon,'Color','k')
% ylim([15 45])
ylim([0 90])
% xlim([-180 -120])
xlim([-280 -100])
% caxis([3 7])
xlabel('Longitude (ºE)')
ylabel('Latitude (ºN)')
cb1 = colorbar;

subplot(2,1,2)
mat2_xy(mat2_xy < thr) = NaN;
contourf(lon_bin, lat_bin, log( mat2_xy ), 100, 'linestyle', 'none')
geoshow(coastlat,coastlon,'Color','k')
% ylim([15 45])
ylim([0 90])
% xlim([-180 -120])
xlim([-280 -100])
% caxis([3 7])
xlabel('Longitude (ºE)')
ylabel('Latitude (ºN)')
cb2 = colorbar;


f1 = figure('Position', [1 1 800 400]);
contourf(lon_bin, lat_bin, ( squeeze( mat2_xy-mat1_xy ) ), 100, 'linestyle', 'none')
geoshow(coastlat,coastlon,'Color','k')
ylim([15 45])
xlim([-180 -120])
xlabel('Longitude (ºE)')
ylabel('Latitude (ºN)')
% caxis([-250 250])
colorbar;
colormap('redblue')
% exportgraphics(f1, 'hysplit_latlon_diff.eps')


%% Time-latitude


close all
f1 = figure('Position', [1 1 1200 500]);
subplot(1,2,1)
mat1_y(mat1_y < thr) = NaN;
contourf(lat_bin(2:end), (0:Nh-1)/24, (mat1_y)', 20)
view(90,90)
% xlim([300 1000])
% xlabel('Pressure (hPa)')
ylabel('Time (Day)')
cb = colorbar;
title(cb, '(#)')

subplot(1,2,2)
mat2_y(mat2_y < thr) = NaN;
contourf(lat_bin(2:end), (0:Nh-1)/24, (mat2_y)', 20)
view(90,90)
% xlim([300 1000])
% xlabel('Pressure (hPa)')
ylabel('Time (Day)')
cb = colorbar;
title(cb, '(#)')


f2 = figure('Position', [1 1 700 500]);
contourf(lat_bin(2:end), (0:Nh-1)/24, (mat1_y-mat2_y)', -200:2:200, 'linestyle', 'none')
view(90,90)
% xlim([300 1000])
% xlabel('Pressure (hPa)')
ylabel('Time (Day)')
colormap('redblue')
caxis([-150 150])
cb = colorbar;
title(cb, '(#)')
% exportgraphics(f2, 'hysplit_prestime_diff.eps')


%% Time - Latitude CDF

cmp = lines(7);

d1 =  25;
d3 =  73;
d5 = 121;

close all
f3 =  figure('Position', [1 1 700 500]);
hold on

% Day 1
mat = cumsum( mat1_y(:,d1))./sum( mat1_y(:,1));
plot(lat_bin(2:end), mat, 'color', cmp(1,:), 'linestyle', ':')

mat = cumsum( mat2_y(:,d1))./sum( mat2_y(:,1));
plot(lat_bin(2:end), mat, 'color', cmp(2,:), 'linestyle', ':')

% Day 2
mat = cumsum( mat1_y(:,d3))./sum( mat1_y(:,1));
plot(lat_bin(2:end), mat, 'color', cmp(1,:), 'linestyle', '--')

mat = cumsum( mat2_y(:,d3))./sum( mat2_y(:,1));
plot(lat_bin(2:end), mat, 'color', cmp(2,:), 'linestyle', '--')

% Day 3
mat = cumsum( mat1_y(:,d5))./sum( mat1_y(:,1));
[~,pos] = min(abs(mat-0.75));
lat_w   = lat_bin(pos);
plot(lat_bin(2:end), mat, 'color', cmp(1,:))

mat = cumsum( mat2_y(:,d5))./sum( mat2_y(:,1));
[~,pos] = min(abs(mat-0.75));
lat_d   = lat_bin(pos);
plot(lat_bin(2:end), mat, 'color', cmp(2,:))

grid on
box on

set(gca, 'ytick', 0:0.25:1)

ylim([0 1])
xlim([15 60])

%% Time - Longitude CDF

cmp = lines(7);

d1 =  25;
d3 =  73;
d5 = 121;

% close all
f3 =  figure('Position', [1 1 700 500]);
hold on

% Day 1
mat = cumsum( mat1_x(:,d1))./sum( mat1_x(:,1));
plot(lon_bin(2:end), mat, 'color', cmp(1,:), 'linestyle', ':')

mat = cumsum( mat2_x(:,d1))./sum( mat2_x(:,1));
plot(lon_bin(2:end), mat, 'color', cmp(2,:), 'linestyle', ':')

% Day 2
mat = cumsum( mat1_x(:,d3))./sum( mat1_x(:,1));
plot(lon_bin(2:end), mat, 'color', cmp(1,:), 'linestyle', '--')

mat = cumsum( mat2_x(:,d3))./sum( mat2_x(:,1));
plot(lon_bin(2:end), mat, 'color', cmp(2,:), 'linestyle', '--')

% Day 3
mat = cumsum( mat1_x(:,d5))./sum( mat1_x(:,1));
[~,pos] = min(abs(mat-0.25));
lon_w   = lon_bin(pos);
plot(lon_bin(2:end), mat, 'color', cmp(1,:))

mat = cumsum( mat2_x(:,d5))./sum( mat2_x(:,1));
[~,pos] = min(abs(mat-0.25));
lon_d   = lon_bin(pos);
plot(lon_bin(2:end), mat, 'color', cmp(2,:))

grid on
box on

set(gca, 'ytick', 0:0.25:1)

ylim([0 1])
xlim([-250 -120])



%% Special plot #1

close all

cmp = lines(7);

d1 =  25;
d3 =  73;
d5 = 121;


f2 = figure('Position', [1 1 800 700]);

sb1 = subplot(2,1,1);
contourf(lon_bin, lat_bin, ( squeeze( mat2_xy-mat1_xy ) ), 100, 'linestyle', 'none')
geoshow(coastlat,coastlon,'Color','k')
ylim([15 45])
xlim([-180 -120])
xlabel('Longitude (ºE)')
ylabel('Latitude (ºN)')
caxis([-2000 2000])
colorbar;
colormap('redblue')

text(-179,43,'a)', 'fontsize', 19)

sb2 = subplot(2,1,2);

hold on

% Day 1
mat = cumsum( mat1_y(:,d1))./sum( mat1_y(:,1));
plot(lat_bin(2:end), mat, 'color', cmp(1,:), 'linestyle', ':')

mat = cumsum( mat2_y(:,d1))./sum( mat2_y(:,1));
plot(lat_bin(2:end), mat, 'color', cmp(2,:), 'linestyle', ':')

% Day 2
mat = cumsum( mat1_y(:,d3))./sum( mat1_y(:,1));
plot(lat_bin(2:end), mat, 'color', cmp(1,:), 'linestyle', '--')

mat = cumsum( mat2_y(:,d3))./sum( mat2_y(:,1));
plot(lat_bin(2:end), mat, 'color', cmp(2,:), 'linestyle', '--')

% Day 3
mat = cumsum( mat1_y(:,d5))./sum( mat1_y(:,1));
plot(lat_bin(2:end), mat, 'color', cmp(1,:))

mat = cumsum( mat2_y(:,d5))./sum( mat2_y(:,1));
plot(lat_bin(2:end), mat, 'color', cmp(2,:))

grid on
box on

set(gca, 'ytick', 0:0.25:1)

ylim([0 1])
xlim([15 60])

xlabel('Latitude (ºN)')
ylabel('CDF')

legend('WET - day 1', 'DRY - day 1', 'WET - day 3',...
    'DRY - day 3', 'WET - day 5', 'DRY - day 5', 'location', 'southeast')


text(15.5,0.93,'b)', 'fontsize', 19)



sb1.Position = [0.0100    0.5438    0.8750    0.3975];
sb2.Position = [0.1000    0.0800    0.7850    0.35];

% exportgraphics(f2, '../plots/hysplit_lat_differences.eps')



%% Special plot #2

close all

cmp = lines(7);

d1 =  25;
d3 =  73;
d5 = 121;


f2 = figure('Position', [1 1 1000 700]);

sb1 = subplot(3,1,1);
contourf(lon_bin, lat_bin, ( squeeze( mat2_xy-mat1_xy ) ), 100, 'linestyle', 'none')
geoshow(coastlat,coastlon,'Color','k')
ylim([15 45])
xlim([-180 -120])
xlabel('Longitude (ºE)')
ylabel('Latitude (ºN)')
caxis([-2000 2000])
colorbar;
colormap('redblue')

text(-179,43,'a)', 'fontsize', 19)


sb2 = subplot(3,1,2);

hold on

% Day 1
mat = cumsum( mat1_x(:,d1))./sum( mat1_x(:,1));
plot(lon_bin(2:end), 1-mat, 'color', cmp(1,:), 'linestyle', ':')

mat = cumsum( mat2_x(:,d1))./sum( mat2_x(:,1));
plot(lon_bin(2:end), 1-mat, 'color', cmp(2,:), 'linestyle', ':')

% Day 2
mat = cumsum( mat1_x(:,d3))./sum( mat1_x(:,1));
plot(lon_bin(2:end), 1-mat, 'color', cmp(1,:), 'linestyle', '--')

mat = cumsum( mat2_x(:,d3))./sum( mat2_x(:,1));
plot(lon_bin(2:end), 1-mat, 'color', cmp(2,:), 'linestyle', '--')

% Day 3
mat = cumsum( mat1_x(:,d5))./sum( mat1_x(:,1));
plot(lon_bin(2:end), 1-mat, 'color', cmp(1,:))

mat = cumsum( mat2_x(:,d5))./sum( mat2_x(:,1));
plot(lon_bin(2:end), 1-mat, 'color', cmp(2,:))

grid on
box on


ylim([0 1])
xlim([-220 -100])

set(gca, 'ytick', 0:0.25:1)
set(gca, 'xtick', -220:40:0)
set(gca, 'xticklabel', [140 160 -180 -160 -140 -120 -100])

xlabel('Longitude (ºE)')
ylabel('CDF')

% text(-218,0.93,'b)', 'fontsize', 19)

text(-109,0.93,'b)', 'fontsize', 19)



sb3 = subplot(3,1,3);

hold on

% Day 1
mat = cumsum( mat1_y(:,d1))./sum( mat1_y(:,1));
plot(lat_bin(2:end), mat, 'color', cmp(1,:), 'linestyle', ':')

mat = cumsum( mat2_y(:,d1))./sum( mat2_y(:,1));
plot(lat_bin(2:end), mat, 'color', cmp(2,:), 'linestyle', ':')

% Day 2
mat = cumsum( mat1_y(:,d3))./sum( mat1_y(:,1));
plot(lat_bin(2:end), mat, 'color', cmp(1,:), 'linestyle', '--')

mat = cumsum( mat2_y(:,d3))./sum( mat2_y(:,1));
plot(lat_bin(2:end), mat, 'color', cmp(2,:), 'linestyle', '--')

% Day 3
mat = cumsum( mat1_y(:,d5))./sum( mat1_y(:,1));
plot(lat_bin(2:end), mat, 'color', cmp(1,:))

mat = cumsum( mat2_y(:,d5))./sum( mat2_y(:,1));
plot(lat_bin(2:end), mat, 'color', cmp(2,:))

grid on
box on

set(gca, 'ytick', 0:0.25:1)
set(gca, 'xtick', 0:15:90)
ylim([0 1])
xlim([15 60])
legend('WET - day 1', 'DRY - day 1', 'WET - day 3',...
    'DRY - day 3', 'WET - day 5', 'DRY - day 5', 'location', 'southeast')

xlabel('Latitude (ºN)')
ylabel('CDF')


text(16,0.93,'c)', 'fontsize', 19)


sb1.Position = [0.0500    0.5438    0.8750    0.3975];
sb2.Position = [0.1000    0.0800    0.375    0.35];
sb3.Position = [0.5700    0.0800    0.375    0.35];

exportgraphics(f2, '../plots/hysplit_lat_lon_differences.eps')



%% Time-pressure


figure('Position', [1 1 1200 500]);
subplot(1,2,1)
mat1_p(mat1_p < thr) = NaN;
contourf(p_bin(2:end), (0:Nh-1)/24, log(mat1_p)', 0:.5:10)
view(90,90)
xlim([300 1000])
xlabel('Pressure (hPa)')
ylabel('Time (Day)')
cb = colorbar;
title(cb, '(#)')

subplot(1,2,2)
mat2_p(mat2_p < thr) = NaN;
contourf(p_bin(2:end), (0:Nh-1)/24, log(mat2_p)', 0:.5:10)
view(90,90)
xlim([300 1000])
xlabel('Pressure (hPa)')
ylabel('Time (Day)')
cb = colorbar;
title(cb, '(#)')


f2 = figure('Position', [1 1 700 500]);
contourf(p_bin(2:end), (0:Nh-1)/24, (mat1_p-mat2_p)', -200:2:200, 'linestyle', 'none')
view(90,90)
xlim([200 1000])
xlabel('Pressure (hPa)')
ylabel('Time (Day)')
colormap('redblue')
caxis([-150 150])
cb = colorbar;
title(cb, '(#)')
% exportgraphics(f2, 'hysplit_prestime_diff.eps')



cmp = lines(7);

d1 =  25;
d3 =  73;
d5 = 121;

% close all
f3 =  figure('Position', [1 1 700 500]);
hold on

% Day 1
mat = cumsum( mat1_p(:,d1))./sum( mat1_p(:,1));
plot(p_bin(2:end), mat, 'color', cmp(1,:), 'linestyle', ':')

mat = cumsum( mat2_p(:,d1))./sum( mat2_p(:,1));
plot(p_bin(2:end), mat, 'color', cmp(2,:), 'linestyle', ':')

% Day 2
mat = cumsum( mat1_p(:,d3))./sum( mat1_p(:,1));
plot(p_bin(2:end), mat, 'color', cmp(1,:), 'linestyle', '--')

mat = cumsum( mat2_p(:,d3))./sum( mat2_p(:,1));
plot(p_bin(2:end), mat, 'color', cmp(2,:), 'linestyle', '--')

% Day 3
mat = cumsum( mat1_p(:,d5))./sum( mat1_p(:,1));
plot(p_bin(2:end), mat, 'color', cmp(1,:))

mat = cumsum( mat2_p(:,d5))./sum( mat2_p(:,1));
plot(p_bin(2:end), mat, 'color', cmp(2,:))

grid on
box on

set(gca, 'ytick', 0:0.25:1)

ylim([0 1])
% xlim([15 45])





