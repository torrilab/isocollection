set(0, 'defaultaxesfontsize',19)
set(0, 'defaultlinelinewidth', 2)


clear all
close all

root = '../data/';

mat  = xlsread([root,'brian_data_new.xlsx']);

i0   = 7;

mat(1:i0,:) = [];

Nd   = size(mat,1);

d0   = datenum('01-Jan-1900');
d1   = datenum('01-Jan-2019');


day  = mat(:,3)+d0-d1-1;

mn   = month(datetime(datestr(day)));

ddt  = mat(:,4);
rml  = mat(:,5);
rmm  = mat(:,7);
mld  = mat(:,6);
mmd  = rmm./ddt;


d18o = mat(:,11);
dD   = mat(:,13);
dexc = dD - 8*d18o;

d18O = d18o;
save([root,'brian_data.mat'], 'Nd', 'day', 'rml', 'rmm', 'mmd', 'd18O', 'dD', 'dexc')

%% Separate data in wet and dry season


vec_l = ~isnan(dD);

mn1     = mn(vec_l);
dD1     = dD(vec_l);
rmm1    = rmm(vec_l);
rml1    = rml(vec_l);
mld1    = mld(vec_l);
mmd1    = mmd(vec_l);
d18o1   = d18o(vec_l);
dexc1   = dexc(vec_l);
% Nd1     = Nd - sum(vec_l);
Nd1     = length(mn1);


vec_d = mn1 > 4 & mn1 <= 10;
vec_w = mn1 <= 4 | mn1 > 10;

mn_d    = mn1(vec_d);
dD_d    = dD1(vec_d);
rml_d   = rml1(vec_d);
mld_d   = mld1(vec_d);
mmd_d   = mmd1(vec_d);
d18o_d  = d18o1(vec_d);
dexc_d  = dexc1(vec_d);

Nd_d    = length(mn_d);

mn_w    = mn1(vec_w);
dD_w    = dD1(vec_w);
rml_w   = rml1(vec_w);
mld_w   = mld1(vec_w);
mmd_w   = mmd1(vec_w);
d18o_w  = d18o1(vec_w);
dexc_w  = dexc1(vec_w);

Nd_w    = length(mn_w);

%% Rainfall


close all
figure('Position', [1 1 1230 400])

yyaxis left
vec = mmd;
vec(vec == 0) = NaN;
plot(day,vec, '-o')
% ylim([0 40])
ylabel('Rainfall rate (mm day^{-1})')

yyaxis right
hold on
vec = rmm;
vec(vec == 0) = NaN;
plot(day, vec, '-o')
hold off
% ylim([0 80])
ylabel('Collected rainfall (mm)')


xlim([190 920])
set(gca, 'xtick', [244 335 426 518 610 701 792 884])
set(gca, 'xticklabel', [{'Sep19', 'Dec19', 'Mar20', 'Jun20', 'Sep20', 'Dec20', 'Mar21', 'Jun20'}])
grid on
box on



%% Isotopic composition + d-excess


close all
figure('Position', [1 1 1230 700])

sb1 = subplot(2,1,1);

yyaxis right
hold on
plot(day, d18o, '-o')
% plot(day, movmean(d18o,4, 'omitnan'))
hold off
ylim([-10 2])
% ylim([-2 2])
ylabel(['\delta^{18}O (',char(8240),')'])

yyaxis left
plot(day,dD, '-o')
ylim([-80 40])
% ylim([-20 20])
ylabel(['\deltaD (',char(8240),')'])

xlim([190 920])
set(gca, 'xtick', [244 335 426 518 610 701 792 884])
set(gca, 'xticklabel', [{'Sep19', 'Dec19', 'Mar20', 'Jun20', 'Sep20', 'Dec20', 'Mar21', 'Jun20'}])
grid on
box on

sb1.Position = [0.0600    0.5838    0.8750    0.4];
text( 193, 34,'a)', 'fontsize', 19) 


sb2 = subplot(2,1,2);
plot(day, dexc, 'k-o')
ylabel(['d-excess (',char(8240),')'])

ylim([6 22])
set(gca, 'ytick', 6:4:22)
xlim([190 920])
set(gca, 'xtick', [244 335 426 518 610 701 792 884])
set(gca, 'xticklabel', [{'Sep19', 'Dec19', 'Mar20', 'Jun20', 'Sep20', 'Dec20', 'Mar21', 'Jun20'}])
grid on
box on

sb2.Position = [0.0600    0.0800    0.8750    0.4];
text( 193, 21.3,'b)', 'fontsize', 19) 


%%


h1 = histogram(dD, -100:1:10, 'normalization', 'cumcount');
cdat_dD   = h1.Values;
cbin_dD   = h1.BinEdges;


h2 = histogram(d18o, -10:.1:2, 'normalization', 'cumcount');
cdat_d18o = h2.Values;
cbin_d18o = h2.BinEdges;

close all


figure('Position', [1 1 700 500])
% subplot(1,2,1)

yyaxis left
plot(cdat_dD/Nd1, cbin_dD(1:end-1))
grid on

% subplot(1,2,2)
yyaxis right
plot(cdat_d18o/Nd1, cbin_d18o(1:end-1))
grid on

% xlim([-70 10])

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

xlim([-70 10])
ylim([0 0.12])
set(gca, 'ytick', 0:0.04:0.12)
box on
grid on
ylabel('PDF')
xlabel(['\deltaD (',char(8240),')'])
legend('WET', 'DRY', 'location', 'northeast')
sb1.Position = [0.0600    0.5838    0.90    0.3412];

text( -69.5, 0.11,'a)', 'fontsize', 19) 



sb2 = subplot(2,1,2);
x = c18o_bin(1:end-1)';
y = [c18o_wet_val; c18o_dry_val]';
bar(x, y, 'grouped')

xlim([-10 0])
ylim([0 1])
set(gca, 'ytick', 0:0.2:1)
box on
grid on
ylabel('PDF')
xlabel(['\delta^{18}O (',char(8240),')'])
sb2.Position = [0.0600    0.1100    0.90    0.3412];

text( -9.93, 0.92,'b)', 'fontsize', 19) 

%% 3D histogram by season and rainfall rate

cD_bin   = -100:2:10;

mmd_bin  = 0:.2:10;


h3D_dry = hist3([dD_d, mmd_d], {cD_bin,mmd_bin})/Nd_d;
h3D_wet = hist3([dD_w, mmd_w], {cD_bin,mmd_bin})/Nd_w;
mat = h3D_wet - h3D_dry;
mat(mat == 0) = NaN;
pcolor(mmd_bin,cD_bin,mat)
caxis([-.1 0.1])
% colormap('parula')
colormap('redblue')

ylim([-40 10])
xlim([0 10])
%%
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






%% Histograms of abundance and rainfall
% I'm not sure about this part because data from Brian has been 
% collected in unequal intervals, so we're weighing longer time collections
% more than the shorter ones


d_iso    = dD1;
r_iso    = rmm1;

bin_iso  = -100:.01:10;

h1       = histogram(d_iso, bin_iso, 'normalization', 'cumcount');
cdat_dD  = h1.Values/Nd1;
cbin_dD  = h1.BinEdges;

cbin_rml = 0:0.001:1;
Nrml     = length(cbin_rml);
cdat_rml = zeros(1,Nrml);
Nnorm    = nansum(r_iso);

for n = 1:Nrml
    
    i = cbin_rml(n);
    
    [~,pos] = min(abs( cdat_dD-i ));
    dD_thr  = cbin_dD(pos);
    idx     = find(d_iso < dD_thr);
    
    cdat_rml(n) = nansum( r_iso(idx) );
        
end


close all
figure('Position', [1 1 700 600])

yyaxis left
plot(cdat_dD, cbin_dD(1:end-1))
ylim([-70 10])
grid on
ylabel(['\deltaD (',char(8240),')'])

yyaxis right
plot(cdat_rml/Nnorm, cbin_rml)
line([0 1],[0 1], 'color', 'k')
set(gca, 'ytick', 0:0.125:1)
ylabel('CDF (prec.)')

set(gca, 'xtick', 0:0.25:1)
xlabel('CDF (\deltaD)')


%% Amount Effect

vin  = mmd1;
iso1 = dD1;
iso2 = d18o1;
Nx   = Nd1;

bin = 0:1:100;
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
ylim([-60 10])
ylabel(['\deltaD (',char(8240),')'])

yyaxis right
scatter(bin(1:end-1),amt(2,:), 100, 'fill', 'markeredgecolor', 'k')
ylabel(['\delta^{18}O (',char(8240),')'])
ylim([-8 -1])

xlim([0 40])
grid on
box on

xlabel('Precipitation rate (mm day^{-1})')






%% Local Meteoric Water Line


% Linear regression of all the data
mn1     = mn(~isnan(dD));
dD1     = dD(~isnan(dD));
d18o1   = d18o(~isnan(d18o));
mmd1    = mmd(~isnan(d18o));
mdl_tot = fitlm(d18o1,dD1);


% Linear regression of dry season data
mn2     = mn1(mn1 >= 4 & mn1 <= 10);
dD2     = dD1(mn1 >= 4 & mn1 <= 10);
d18o2   = d18o1(mn1 >= 4 & mn1 <= 10);
mmd2    = mmd(mn1 >= 4 & mn1 <= 10);

mdl_dry = fitlm(d18o2,dD2);


% Linear regression of wet season data
mn3     = mn1(mn1 < 4 | mn1 > 10);
dD3     = dD1(mn1 < 4 | mn1 > 10);
d18o3   = d18o1(mn1 < 4 | mn1 > 10);
mmd3    = mmd(mn1 < 4 | mn1 > 10);
mdl_wet = fitlm(d18o3,dD3);



% Plot regressions
% cmp     = cat(1, parula(6), flip(parula(6)));

% v_wet = [0.2422    0.1504    0.6603];
% v_dry = [0.9769    0.9839    0.0805];

vec   = lines(2);
v_wet = vec(1,:);
v_dry = vec(2,:);

cmp1 = repmat(v_wet,[4,1]);
cmp2 = repmat(v_dry,[6,1]);
cmp3 = repmat(v_wet,[2,1]);

cmp = cat(1,cmp1,cmp2,cmp3);

close all
figure('Position', [1 1 600 600])
c = scatter(d18o1, dD1, 150, cmp(mn1,:), 'filled');
c.MarkerEdgeColor = 'k';

ylabel(['\deltaD (',char(8240),')'])
xlabel(['\delta^{18}O (',char(8240),')'])

grid on
box on

xlim([-10 2])
ylim([-80 20]) 

% colormap(cmp)
% cb = colorbar;
% cb.Ticks = (0.5:11.5)/12;
% cb.TickLabels = 1:12;

% title('Local Meteoric Water Line')




%% Amount effect by season


close all
figure('Position', [1 1 1200 600])
subplot(1,2,1)
hold on

c1 = scatter(mmd2, d18o2, 150, [0.8500    0.3250    0.0980], 'filled');
c1.MarkerEdgeColor = 'k';

c2 = scatter(mmd3, d18o3, 150, [0    0.4470    0.7410], 'filled');
c2.MarkerEdgeColor = 'k';

xlim([0 16])
set(gca, 'xtick', 0:4:16)
ylabel(['\delta^{18}O (',char(8240),')'])
xlabel('Rain rate (mm day^{-1})')

grid on
box on


subplot(1,2,2)
hold on

c1 = scatter(mmd2, dD2, 150, [0.8500    0.3250    0.0980], 'filled');
c1.MarkerEdgeColor = 'k';

c2 = scatter(mmd3, dD3, 150, [0    0.4470    0.7410], 'filled');
c2.MarkerEdgeColor = 'k';

xlim([0 16])
set(gca, 'xtick', 0:4:16)
ylabel(['\deltaD (',char(8240),')'])
xlabel('Rain rate (mm day^{-1})')

grid on
box on


