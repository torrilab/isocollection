% Needs nino_indices.m

% clear all

set(0, 'DefaultAxesFontSize', 19)
set(0, 'DefaultLineLineWidth', 2)

roots = '../data/';

load([roots,'pacific_indices.mat'])
% Download the data

% MAT   = readtable([roots,'gnip_hilo.xlsx']);
MAT   = xlsread([roots,'gnip_hilo.xlsx']);


yr_i  = MAT(:,12);
yr_f  = MAT(:,13);
d18o  = MAT(:,14);
dD    = MAT(:,16);
dT    = MAT(:,18);
prec  = MAT(:,21);

dexc = dD - 8*d18o;


% yr0   = year(yr_i(1));
Ndata  = size(MAT,1);
d_orig = yr_i(1);
d0     = datenum('01-01-1960');
yr_i   = yr_i-d_orig+d0;
yr_f   = yr_f-d_orig+d0;



mat   = mat_indices(:,3);
enso  = mat(:);

mat   = mat_indices(:,4);
pdo   = mat(:);

mat   = mat_indices(:,5);
pmm   = mat(:);



%% Find months and years of data

tmp1 = ( datetime(datestr(yr_i)) );
tmp2 = ( datetime(datestr(yr_f)) );

yr_i = year(tmp1);
mn_i = month(tmp1);

yr_f = year(tmp2);
mn_f = month(tmp2);

dm   = eomday(yr_i,mn_i);
prate  = prec./dm;


%% Divide data into wet and dry season

vec_l   = ~isnan(dD) & ~isnan(d18o);

mn1     = mn_i(vec_l);
yr1     = yr_i(vec_l);
dm1     = dm(vec_l);
dD1     = dD(vec_l);
prec1   = prec(vec_l);
prate1  = prec1./dm1;
d18o1   = d18o(vec_l);
EL3     = enso(vec_l);
PDO     = pdo(vec_l);
PMM     = pmm(vec_l);
Nd1     = length(mn1);
dexc1   = dD1-8*d18o1;


vec_d   = mn1 > 4 & mn1 <= 10;
vec_w   = mn1 <= 4 | mn1 > 10;


mn_d    = mn1(vec_d);
yr_d    = yr1(vec_d);
dD_d    = dD1(vec_d);
prec_d  = prec1(vec_d);
prate_d = prate1(vec_d);
d18o_d  = d18o1(vec_d);
EL3_d   = enso(vec_d);
PDO_d   = pdo(vec_d);
PMM_d   = pmm(vec_d);
Nd_d    = length(mn_d);
dexc_d  = dD_d-8*d18o_d;


mn_w    = mn1(vec_w);
yr_w    = yr1(vec_w);
dD_w    = dD1(vec_w);
prec_w  = prec1(vec_w);
prate_w = prate1(vec_w);
d18o_w  = d18o1(vec_w);
EL3_w   = enso(vec_w);
PDO_w   = pdo(vec_w);
PMM_w   = pmm(vec_w);
Nd_w    = length(mn_w);
dexc_w  = dD_w-8*d18o_w;


%% Rainfall rates

close all
figure('Position', [1 1 1230 400])




yyaxis right

plot(prate, '-o')




yyaxis left
plot(dD, '-o')
ylim([-30 10])
set(gca, 'ytick', -30:5:10)
ylabel(['\delta^{2}H (',char(8240),')'])

xlim([24 120])
set(gca, 'xtick', 0:12:144)
set(gca, 'xticklabel', [{'1960', '1961', '1962', '1963', '1964', '1965', ...
     '1966', '1967', '1968', '1969', '1970'}])
grid on
box on



% set(gca, 'ytick', -7:1:1)
% ylabel(['\delta^{18}O (',char(8240),')'])
% 
% 
% % xlim([24 120])
% set(gca, 'xtick', 0:12:144)
% set(gca, 'xticklabel', [{'1960', '1961', '1962', '1963', '1964', '1965', ...
%      '1966', '1967', '1968', '1969', '1970'}])
% grid on
% box on




%% Isotopic composition + d-excess

close all
figure('Position', [1 1 1230 700])

sb1 = subplot(2,1,1);

yyaxis right
plot(d18o, '-o')
ylim([-7 1])
% set(gca, 'ytick', -7:1:1)
ylabel(['\delta^{18}O (',char(8240),')'])


yyaxis left
plot(dD, '-o')
ylim([-30 10])
set(gca, 'ytick', -30:5:10)
ylabel(['\delta^{2}H (',char(8240),')'])

xlim([24 120])
set(gca, 'xtick', 0:12:144)
set(gca, 'xticklabel', [{'1960', '1961', '1962', '1963', '1964', '1965', ...
     '1966', '1967', '1968', '1969', '1970'}])
grid on
box on

sb1.Position = [0.0600    0.5838    0.8750    0.4];
text( 24.5, 8,  'a)', 'fontsize', 19) 


sb2 = subplot(2,1,2);
plot(dexc, 'k-o')
ylabel(['d-excess (',char(8240),')'])

xlim([24 120])
set(gca, 'xtick', 0:12:144)
set(gca, 'xticklabel', [{'1960', '1961', '1962', '1963', '1964', '1965', ...
     '1966', '1967', '1968', '1969', '1970'}])
grid on
box on

sb2.Position = [0.0600    0.0800    0.8750    0.4];
text( 24.5, 23,'b)', 'fontsize', 19) 



%% Histograms by season

cD_bin   = -100:4:10;
c18o_bin = -10:.5:1;

% hD_tot   = histogram(dD1,  cD_bin, 'normalization', 'pdf');
% cD_tot_val = hD_tot.Values;

hD_dry   = histogram(dD_d, cD_bin, 'normalization', 'pdf');
cD_dry_val = hD_dry.Values;

hD_wet   = histogram(dD_w, cD_bin, 'normalization', 'pdf');
cD_wet_val = hD_wet.Values;


% hD_tot   = histogram(d18o1,  c18o_bin, 'normalization', 'pdf');
% c18o_tot_val = h18o_tot.Values;

h18o_dry = histogram(d18o_d, c18o_bin, 'normalization', 'pdf');
c18o_dry_val = h18o_dry.Values;

h18o_wet = histogram(d18o_w, c18o_bin, 'normalization', 'pdf');
c18o_wet_val = h18o_wet.Values;



close all
figure('Position', [1 1 1200 600])

sb1 = subplot(2,1,1);
x = cD_bin(1:end-1)';
y = [cD_wet_val; cD_dry_val]';
bar(x, y, 'grouped')

xlim([-40 10])
ylim([0 0.12])
set(gca, 'ytick', 0:0.04:0.12)
box on
grid on
ylabel('PDF')
xlabel(['\delta^{2}H (',char(8240),')'])
legend('WET', 'DRY', 'location', 'northeast')
sb1.Position = [0.0600    0.5838    0.90    0.3412];

text( -69.5, 0.11,'a)', 'fontsize', 19) 



sb2 = subplot(2,1,2);
x = c18o_bin(1:end-1)';
y = [c18o_wet_val; c18o_dry_val]';
bar(x, y, 'grouped')

xlim([-6 1])
ylim([0 1])
set(gca, 'ytick', 0:0.2:1)
box on
grid on
ylabel('PDF')
xlabel(['\delta^{18}O (',char(8240),')'])
sb2.Position = [0.0600    0.1100    0.90    0.3412];

text( -9.93, 0.92,'b)', 'fontsize', 19) 



%% Amount Effect

vin  = prate1;
iso1 = dD1;
iso2 = d18o1;
Nx   = Nd1;
bin  = 0:1:50;

h   = histogram(vin, bin);
vec = h.BinEdges;
Nh  = length(vec)-1;
dv  = h.BinWidth;

vx  = (vin-mod(vin,dv))/dv + 1;
amt = zeros(2,Nh);
num = zeros(1,Nh);

for nd = 1:Nx
    i = vx(nd);
    amt(1,i) = amt(1,i) + iso1(nd);
    amt(2,i) = amt(2,i) + iso2(nd);
    num(i)   = num(i) + 1;
end
amt(1,:) = amt(1,:)./num(:)';
amt(2,:) = amt(2,:)./num(:)';


close all
figure('Position', [1 1 600 500])

yyaxis left
scatter(bin(1:end-1),amt(1,:), 100, 'd', 'fill', 'markeredgecolor', 'k')
ylim([-30 10])
ylabel(['\delta^{2}H (',char(8240),')'])

yyaxis right
scatter(bin(1:end-1),amt(2,:), 100, 'fill', 'markeredgecolor', 'k')
ylabel(['\delta^{18}O (',char(8240),')'])
ylim([-5 -1])

xlim([0 30])
grid on
box on

xlabel('Precipitation rate (mm day^{-1})')


%% Local Meteoric Water Line


% Linear regression of all the data
mdl_tot = fitlm(d18o1,dD1);

% Linear regression of dry season data
mdl_dry = fitlm(d18o_d,dD_d);

% Linear regression of wet season data
mdl_wet = fitlm(d18o_w,dD_w);


% Plot regressions
vec   = lines(2);
v_wet = vec(1,:);
v_dry = vec(2,:);

cmp1 = repmat(v_wet,[4,1]);
cmp2 = repmat(v_dry,[6,1]);
cmp3 = repmat(v_wet,[2,1]);

cmp = cat(1,cmp1,cmp2,cmp3);

close all
figure('Position', [1 1 600 600])
c = scatter(dD1, d18o1, 150, cmp(mn1,:), 'filled');
c.MarkerEdgeColor = 'k';
xlabel(['\delta D (',char(8240),')'])
ylabel(['\delta^{18}O (',char(8240),')'])

grid on
box on

colormap(cmp)
cb = colorbar;
cb.Ticks = (0.5:11.5)/12;
cb.TickLabels = 1:12;

title('Local Meteoric Water Line')


%% LMWL by year

% tmp_w = PDO_w;
% tmp_d = PDO_d;

tmp_w = EL3_w;
tmp_d = EL3_d;

set(0, 'DefaultAxesFontSize', 15)

close all
f = figure('Position', [1 1 1000 800]);
subplot(2,2,1)
hold on
b = scatter(d18o_w, dD_w, 150, (yr_w+mn_w/12), 'filled');
b.MarkerEdgeColor = 'k';
% title('Wet season')
grid on
box on
ylim([-30 10])
xlim([-5 1])
caxis([1962 1970])
ylabel(['\delta^{2}H (',char(8240),')'])
xlabel(['\delta^{18}O (',char(8240),')'])
cb1 = colorbar;
text(-4.85,8,'a)', 'fontsize', 15)
title('Wet Season')

subplot(2,2,2)
hold on
c = scatter(d18o_d, dD_d, 150, (yr_d+mn_d/12), 'filled');
c.MarkerEdgeColor = 'k';
grid on
box on
ylim([-30 10])
xlim([-5 1])
caxis([1962 1970])
ylabel(['\delta^{2}H (',char(8240),')'])
xlabel(['\delta^{18}O (',char(8240),')'])
cb2 = colorbar;
text(-4.85,8,'b)', 'fontsize', 15)
title('Dry Season')

subplot(2,2,3)
hold on
b = scatter(d18o_w, dexc_w, 150, (yr_w+mn_w/12), 'filled');
b.MarkerEdgeColor = 'k';
grid on
box on
xlim([-5 1])
ylim([0 25])
caxis([1962 1970])
xlabel(['\delta^{18}O (',char(8240),')'])
ylabel(['d-excess (',char(8240),')'])
cb3 = colorbar;
text(-4.85,23.5,'c)', 'fontsize', 15)


subplot(2,2,4)
hold on
c = scatter(d18o_d, dexc_d,  150, (yr_d+mn_d/12), 'filled');
c.MarkerEdgeColor = 'k';
grid on
box on
xlim([-5 1])
ylim([0 25])
caxis([1962 1970])
xlabel(['\delta^{18}O (',char(8240),')'])
ylabel(['d-excess (',char(8240),')'])
cb4 = colorbar;
text(-4.85,23.5,'d)', 'fontsize', 15)

exportgraphics(f, '../plots/lmwl_hilo.eps')


%%

set(0, 'DefaultAxesFontSize', 19)

Nyr = length(yr1)/12;

mdl_yr = cell(1,Nyr);

alp = zeros(1,Nyr);
bet = zeros(1,Nyr);


for nyr = 1:Nyr
    
    mi = (nyr-1)*12+1;
    mf = nyr*12;
    
    v1 = d18o1(mi:mf);
    v2 = dD1(mi:mf);
    
    mdl_yr{nyr} = fitlm(v1,v2);
        
    mat = table2array( mdl_yr{nyr}.Coefficients );
    
    alp(nyr) = mat(1,1);
    bet(nyr) = mat(2,1);
    
    disp(nyr)
    
end


%%

v_yr = unique(yr1);
Nyr  = length(v_yr);
mdl_yr = cell(1,Nyr);

alp  = zeros(1,Nyr);
bet  = zeros(1,Nyr);


for nyr = 1:Nyr
    
    
    idx = find(yr1 == v_yr(nyr));
%     idx = find(yr1 > 0);
    
    v1 = d18o1(idx);
    v2 = dD1(idx);
    
    mdl_yr{nyr} = fitlm(v1,v2);
        
    mat = table2array( mdl_yr{nyr}.Coefficients );
    
    alp(nyr) = mat(1,1);
    bet(nyr) = mat(2,1);
    
    disp(nyr)
    
end






%% Group by PDO and NINO3.4

vec_1  = EL3;
vec_2  = PDO;
data1  = d18o1;
data2  = dD1;
data3  = dexc1;

dv1    = 0.5;
min_v1 = -3;
max_v2 = 3;
bin_1  = min_v1:dv1:max_v2;
Nb1    = length(bin_1);

dv2    = 0.5;
min_v2 = -3;
max_v2 = 3;
bin_2  = min_v2:dv2:max_v2;
Nb2    = length(bin_2);


Nd     = length(data1);

h3id1   = zeros(Nb1,Nb2);
h3id2   = zeros(Nb1,Nb2);
h3id3   = zeros(Nb1,Nb2);
num     = zeros(Nb1,Nb2);


for nd = 1:Nd

    b1 = vec_1(nd)-min_v1;
    b1 = floor( b1/dv1 ) + 1;

    b2 = vec_2(nd)-min_v2;
    b2 = floor( b2/dv2 ) + 1;
    
    h3id1(b1,b2) = h3id1(b1,b2) + data1(nd);
    h3id2(b1,b2) = h3id2(b1,b2) + data2(nd);
    h3id3(b1,b2) = h3id3(b1,b2) + data3(nd);
    num(b1,b2)  = num(b1,b2) + 1;
    
end

close all

f = figure('Position', [1 1 1200 500]);

subplot(1,2,1)
pcolor(bin_2, bin_1, h3id1./num)

xlabel('PDO index')
ylabel('Nino 3.4 index')

ylim([-2.5 2.5])
xlim([-2.5 2.5])
caxis([-4 -1])
cb1 = colorbar;
title(cb1, ['\delta^{18}O (',char(8240),')'])
text(-2.43,2.35,'a)', 'fontsize', 19)



subplot(1,2,2)
pcolor(bin_2, bin_1, h3id2./num)

xlabel('PDO index')
ylabel('Nino 3.4 index')

ylim([-2.5 2.5])
xlim([-2.5 2.5])
caxis([-20 2])
% caxis([5 15])
cb2 = colorbar;


title(cb2, ['\delta^{2}H (',char(8240),')'])
text(-2.43,2.35,'b)', 'fontsize', 19)

% exportgraphics(f, '../plots/nino_pdo.eps')


%%


vec_1  = EL3;
vec_2  = PDO;
vec_3  = PMM;

dv1    = 0.5;
min_v1 = -3;
max_v2 = 3;
bin_1  = min_v1:dv1:max_v2;
Nb1    = length(bin_1);

dv2    = 0.5;
min_v2 = -3;
max_v2 = 3;
bin_2  = min_v2:dv2:max_v2;
Nb2    = length(bin_2);

dv3    = 1;
min_v3 = -10;
max_v3 = 10;
bin_3  = min_v3:dv3:max_v3;
Nb3    = length(bin_3);




dat_w  = d18o_w;
dat_d  = d18o_d;

Nd_w   = length(dat_w);
Nd_d   = length(dat_d);

h3id_w = zeros(Nb1,Nb2);
h3id_d = zeros(Nb1,Nb2);
num_w  = zeros(Nb1,Nb2);
num_d  = zeros(Nb1,Nb2);


for nd = 1:Nd_w

    b1 = vec_1(nd)-min_v1;
    b1 = floor( b1/dv1 ) + 1;

    b2 = vec_2(nd)-min_v2;
    b2 = floor( b2/dv2 ) + 1;
    
    h3id_w(b1,b2) = h3id_w(b1,b2) + dat_w(nd);
    num_w(b1,b2)  = num_w(b1,b2) + 1;
    
end


for nd = 1:Nd_d

    b1 = vec_1(nd)-min_v1;
    b1 = floor( b1/dv1 ) + 1;

    b2 = vec_2(nd)-min_v2;
    b2 = floor( b2/dv2 ) + 1;
    
    h3id_d(b1,b2) = h3id_d(b1,b2) + dat_d(nd);
    num_d(b1,b2)  = num_d(b1,b2) + 1;
    
end

close all

f = figure('Position', [1 1 1200 500]);

subplot(1,2,1)
pcolor(bin_2, bin_1, h3id_w./num_w)

xlabel('PDO index')
ylabel('Nino 3.4 index')

ylim([-2.5 2.5])
xlim([-2.5 2.5])
caxis([-4 -1])
cb1 = colorbar;
title(cb1, ['\delta^{18}O (',char(8240),')'])
text(-2.43,2.35,'a)', 'fontsize', 19)
title('Wet Season')

subplot(1,2,2)
pcolor(bin_2, bin_1, h3id_d./num_d)

xlabel('PDO index')
ylabel('Nino 3.4 index')

ylim([-2.5 2.5])
xlim([-2.5 2.5])
caxis([-4 -1])
cb2 = colorbar;

title(cb2, ['\delta^{18}O (',char(8240),')'])
text(-2.43,2.35,'b)', 'fontsize', 19)
title('Dry Season')


exportgraphics(f, '../plots/nino_pdo_seasons.eps')


%%

d18o_avg = zeros(1,12);
dD_avg   = zeros(1,12);
dexc_avg = zeros(1,12);

for n = 1:12
    
    idx = find(mn1 == n);
    
    tmp = sum(d18o1(idx).*prec1(idx))./sum(prec1(idx));
    d18o_avg(n) = tmp;
    
    tmp = sum(dD1(idx).*prec1(idx))./sum(prec1(idx));
    dD_avg(n) = tmp;

    tmp = sum(dexc1(idx).*prec1(idx))./sum(prec1(idx));
    dexc_avg(n) = tmp;

end



%%

close all
figure
subplot(1,2,1)
plot(bin_2, nansum(h3id_w,1)./nansum(num_w,1))
title('PDO')

subplot(1,2,2)
plot(bin_1, sum(h3id_d,2)./sum(num_d,2))
title('NINO 3.4')
