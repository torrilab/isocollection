set(0, 'defaultaxesfontsize',19)
set(0, 'defaultlinelinewidth', 2)

clear all
close all

root     = '../data/';
locs     = {'waikiki', 'hig', 'lyon', 'maunawili', 'brian'};
Nlocs    = length(locs);

Nd_all   = cell(1,Nlocs);
day_all  = cell(1,Nlocs);
yr_all   = cell(1,Nlocs);
mn_all   = cell(1,Nlocs);
rml_all  = cell(1,Nlocs);
rmm_all  = cell(1,Nlocs);
mmd_all  = cell(1,Nlocs);
d18O_all = cell(1,Nlocs);
dD_all   = cell(1,Nlocs);
dexc_all = cell(1,Nlocs);

for ns = 1:Nlocs
    
    load([root,char(locs(ns)),'_data.mat'])
    
    Nd_all{ns}   = Nd;
    day_all{ns}  = day;
    yr_all{ns}   = year(datetime(datestr(day)));
    mn_all{ns}   = month(datetime(datestr(day)));

    rml_all{ns}  = rml;
    rmm_all{ns}  = rmm;
    mmd_all{ns}  = mmd;
    d18O_all{ns} = d18O;
    dD_all{ns}   = dD;
    dexc_all{ns} = dexc;
    
end


clear Nd day rml rmm mmd d18O dD dexc

%%

d0 = datenum('01-January-2019');

d_min = 182;
d_max = 975;

deploy = [277, 271, 206, 278,187];

vec_m = {'January','February','March','April','May',...
    'June','July','August','September','October','November','December'};
vec_y = {'2019', '2020', '2021'};

vec_date = zeros(1,3*12);

for yy = 1:3
for mm = 1:12
    
n = mm+12*(yy-1);
vec_date(n) = datenum(['01-',vec_m{mm},'-',vec_y{yy}])-d0+1;
disp(n)

end
end

vec_date( vec_date < d_min ) = [];
vec_date( vec_date > d_max ) = [];

vec_date_name = [{'','','Sep19','' ,'','Dec19','' ,'','Mar20','' ,'', ...
    'Jun20','','','Sep20','','','Dec20','','','Mar21','','','Jun21', '', '', 'Sep21'}];

%% Rainfall

db1       = .01;
bin1      = 0:db1:100;
Nb1       = length(bin1);
hc1d      = zeros(Nb1-1,Nlocs-1);

cmp       = colormap(lines(Nlocs));
cmp1      = cmp;
cmp1(3,:) = [];


close all
f = figure('Position', [1 1 1230 700]);

sb1 = subplot(2,1,1);
hold on
ns = 0;
for n = 1:Nlocs
    
    vec0 = cell2mat(day_all(n));
    vec1 = cell2mat(mmd_all(n));
    vec1(vec1 < 0) = NaN;
    nc = sum(~isnan(vec1));
    ns = ns+1;
    hc1d(:,ns) = histcounts(vec1, bin1,  'normalization', 'cumcount')/nc;
    if (n == 4 || n == 5)
    plot(vec0,vec1, '-o')
    end
    
end
grid on
box on
xlabel('Time')
ylabel('Rain rate (mm day^{-1})')
xlim([182 975])
% ylim([0 40])
ylim([0 20])
set(gca, 'xtick', vec_date)
set(gca, 'xticklabel', vec_date_name)
text(190, 37, 'a)', 'fontsize', 19)


sb2 = subplot(2,1,2);
hold on
tmp = cat(1,[0 0 0 0 0], hc1d);
plot(bin1,tmp)

xlim([0 30])
set(gca, 'xtick', 0:3:30)
set(gca, 'ytick', 0:0.25:1)
grid on
box on

xlabel('Rain rate (mm day^{-1})')
ylabel('CDF')
legend('Waikiki', 'HIG', 'Lyon', 'Maunawili', 'Kailua', 'location', 'southeast')
text(0.2, .92, 'b)', 'fontsize', 19)

% exportgraphics(f, '../plots/rainfall.eps')

for ns = 1:Nlocs
[~,pos] = min(abs(hc1d(:,ns) - 0.5));
disp(bin1(pos))
end

%% Plot isotopic abundances

cmp = lines(5);

d0 = datenum('01-Jan-2019');
d1 = datenum('01-Oct-2019')-d0+1;
d2 = datenum('01-Apr-2020')-d0+1;
d3 = datenum('01-Oct-2020')-d0+1;
d4 = datenum('01-Apr-2021')-d0+1;


xx = datenum('01-Sept-2019')-d0+1;

close all
f1 = figure('Position', [1 1 1230 600]);

sb1 = subplot(2,1,1);
hold on
for ns = 1:Nlocs
    vec0 = cell2mat(day_all(ns));
    vec1 = cell2mat(d18O_all(ns));
    vec2 = cell2mat(mmd_all(ns));
    idx  = find(vec2 < 0);
    plot(vec0,vec1, '-o', 'color', cmp(ns,:))
end
for ns = 1:Nlocs
    vec0 = cell2mat(day_all(ns));
    vec1 = cell2mat(d18O_all(ns));
    vec2 = cell2mat(mmd_all(ns));
    idx  = find(vec2 < 0);
    scatter(vec0(idx), vec1(idx), 300, 'k', 'x', 'linewidth', 3)
end
grid on
box on
ylabel(['\delta^{18}O (',char(8240),')'])
xlim([vec_date(1) vec_date(end)])
ylim([-10 5])
set(gca, 'ytick', -10:2.5:6)
set(gca, 'xtick', vec_date)
set(gca, 'xticklabel', vec_date_name)
sb1.Position = [0.0600    0.5838    0.8750    0.4];
text( 193, 4.0,'a)', 'fontsize', 19) 


sb2 = subplot(2,1,2);
hold on
for ns = 1:Nlocs
    vec0 = cell2mat(day_all(ns));
    vec1 = cell2mat(dD_all(ns));
    vec2 = cell2mat(mmd_all(ns));
    idx  = find(vec2 < 0);
    plot(vec0,vec1, '-o', 'color', cmp(ns,:))
end

for ns = 1:Nlocs
    vec0 = cell2mat(day_all(ns));
    vec1 = cell2mat(dD_all(ns));
    vec2 = cell2mat(mmd_all(ns));
    idx  = find(vec2 < 0);
    scatter(vec0(idx), vec1(idx), 300, 'k', 'x', 'linewidth', 3)
end
grid on
box on
ylabel(['\delta^{2}H (',char(8240),')'])
xlim([vec_date(1) vec_date(end)])
set(gca, 'xtick', vec_date)
set(gca, 'xticklabel', vec_date_name)

ylim([-64 32])
set(gca, 'ytick', -80:16:32)

sb2.Position = [0.0600    0.0800    0.8750    0.4];
text( 193, 25,'b)', 'fontsize', 19) 
legend('Waikiki', 'HIG', 'Lyon', 'Maunawili', 'Kailua', 'location', 'southeast')
% exportgraphics(f1, '../plots/isoabundances_all.eps')

%% Student t-test

x = cell2mat(dD_all(4));
y = cell2mat(dD_all(5));

[h,p,ci,stats] = ttest2(x,y)



%% Plot deuterium excess

cmp = lines(5);

% close all
f1 = figure('Position', [1 1 1230 350]);
hold on
for ns = 1:Nlocs
    vec0 = cell2mat(day_all(ns));
    vec1 = cell2mat(dexc_all(ns));
    vec2 = cell2mat(mmd_all(ns));
    idx  = find(vec2 < 0);
%     if(ns == 4 || ns == 5)
    p1 = plot(vec0,vec1, '-o', 'color', cmp(ns,:));
%     end
end
% for ns = 1:Nlocs
%     vec0 = cell2mat(day_all(ns));
%     vec1 = cell2mat(dexc_all(ns));
%     vec2 = cell2mat(mmd_all(ns));
%     idx  = find(vec2 < 0);
%     scatter(vec0(idx), vec1(idx), 300, 'k', 'x', 'linewidth', 3)    
% end
xlim([vec_date(1) vec_date(end)])
set(gca, 'xtick', vec_date)
set(gca, 'xticklabel', vec_date_name)
ylim([0 24])
set(gca, 'ytick', 0:4:24)
ylabel(['d-excess (',char(8240),')'])
grid on
box on
legend('Waikiki', 'HIG', 'Lyon', 'Maunawili', 'Kailua', 'location', 'northwest')
% exportgraphics(f1, '../plots/dexc_all.eps')



%% Plot amount effect

cmp = colormap(lines(Nlocs));

db1        = 5;
bin1      = 0:db1:1200;
bin2      = 3*bin1;
Nb1        = length(bin1);
hc1d       = zeros(Nb1-1,Nlocs-1);
hc2d       = zeros(Nb1-1,1);

cmp       = colormap(lines(Nlocs));
cmp1      = cmp;
cmp1(3,:) = [];
cmp2      = cmp(3,:);

ameff1    = zeros(Nb1,Nlocs);
numb1     = zeros(Nb1,Nlocs);
ameff2    = zeros(Nb1,1);
numb2     = zeros(Nb1,1);

close all
figure('Position', [1 1 1200 500])
subplot(1,2,1)
hold on
for ns = 1:Nlocs
    
    vec1 = cell2mat(mmd_all(ns));
    vec1(vec1 < 0) = NaN;
    
    vec2 = cell2mat(d18O_all(ns));
    
    scatter(vec1, vec2, 120, cmp(ns,:), 'fill', 'markerEdgecolor', 'k')
    
    for n = 1:length(vec2)
    nr = floor(vec1(n)/db1)+1;
    nr(isnan(nr)) = Nb1;
    if(~isnan(vec2(n)))
    ameff1(nr,ns) = ameff1(nr,ns) + vec2(n);
    numb1(nr,ns)  = numb1(nr,ns)  + 1;
    end
%     disp([nr ns vec3(n)])
    end
        
    
end
xlim([0 60])
set(gca, 'xtick', 0:10:60)
% ylabel(['\delta^{2}H (',char(8240),')'])
ylabel(['\delta^{18}O (',char(8240),')'])
xlabel('Rain rate (mm day^{-1})')
grid on
box on
legend('Waikiki', 'HIG', 'Lyon', 'Maunawili', 'Kailua', 'location', 'southeast')
% text(280, 16, 'a)', 'fontsize', 19)
text(-40, 4, 'a)', 'fontsize', 19)


subplot(1,2,2)
ameff1 = ameff1./numb1;
bar(bin1+db1/2,ameff1)
grid on
xlim([0 60])
legend('Waikiki', 'HIG', 'Lyon', 'Maunawili', 'Kailua', 'location', 'southwest')
% ylabel(['\delta^{2}H (',char(8240),')'])
ylabel(['\delta^{18}O (',char(8240),')'])

xlabel('Rain rate (mm day^{-1})')
% text(1, 7, 'b)', 'fontsize', 19)
text(-9, 0, 'b)', 'fontsize', 19)


%% Plot amount effect by percentiles


db1       = 0.2;
bin1      = 0:db1:1;

db2       = 0.01;
bin2      = 0:db2:100;

Nb1       = length(bin1);
hcpd      = zeros(Nb1-1,Nlocs-1);

cmp       = colormap(lines(Nlocs));

ameff1    = zeros(Nb1,Nlocs);
numb1     = zeros(Nb1,Nlocs);
ameff2    = zeros(Nb1,1);
numb2     = zeros(Nb1,1);

close all
f1 = figure('Position', [1 1 1200 500]);

for ns = 1:Nlocs
    
    
    vec1  = cell2mat(mmd_all(ns));
    vec1(vec1 < 0) = NaN;
    
    h_tmp = histcounts(vec1,bin2, 'normalization', 'cdf');
    v_tmp = floor( (vec1-min(bin2))/db2 )+1;  
    v1    = ~isnan(v_tmp);
    v_pct = h_tmp(v_tmp(v1));
    
    vec2 = cell2mat(d18O_all(ns));
    vec2 = vec2(v1);
    
    
    
    for n = 1:length(vec2)
        
    nr = floor(v_pct(n)/db1)+1;
    nr(isnan(nr)) = Nb1;
    
    if(~isnan(vec2(n)))
    ameff1(nr,ns) = ameff1(nr,ns) + vec2(n);
    numb1(nr,ns)  = numb1(nr,ns)  + 1;
    end
    
    end
    
    
    
    mn  = mn_all{ns};
    v_d = mn > 4 & mn <= 10;
    v_w = mn <= 4 | mn > 10;
    
    v_d = v_d(v1);
    v_w = v_w(v1);
    
    
    v_pctw = v_pct(v_w);
    v_pctd = v_pct(v_d);
    
    vec2w = vec2(v_w);
    vec2d = vec2(v_d);
        
    sb1 = subplot(1,3,1);
    hold on
    scatter(v_pct, vec2, 120, cmp(ns,:), 'fill', 'markerEdgecolor', 'k')
    if(ns == Nlocs)
    xlim([0 1])
    ylim([-10 4])
    set(gca, 'xtick', 0:0.25:1)
    ylabel(['\delta^{18}O (',char(8240),')'])
    xlabel('Rain rate percentile')
    grid on
    box on
    text(0.03, 3.5, 'a)', 'fontsize', 19)
    end
    title('Full Collection Period')
    
    sb2 = subplot(1,3,2);
    hold on
    scatter(v_pctw, vec2w, 120, cmp(ns,:), 'fill', 'markerEdgecolor', 'k')
    if(ns == Nlocs)
    xlim([0 1])
    ylim([-10 4])
    set(gca, 'xtick', 0:0.25:1)
    ylabel(['\delta^{18}O (',char(8240),')'])
    xlabel('Rain rate percentile')

    grid on
    box on
    text(0.03, 3.5, 'b)', 'fontsize', 19)
    end
    title('Wet Season')
    
    sb3 = subplot(1,3,3);
    hold on
    scatter(v_pctd, vec2d, 120, cmp(ns,:), 'fill', 'markerEdgecolor', 'k')
    if(ns == Nlocs)
    xlim([0 1])
    ylim([-10 4])
    set(gca, 'xtick', 0:0.25:1)
    ylabel(['\delta^{18}O (',char(8240),')'])
    xlabel('Rain rate percentile')
    grid on
    box on
    text(0.03, 3.5, 'c)', 'fontsize', 19)
    legend('Waikiki', 'HIG', 'Lyon', 'Maunawili', 'Kailua', 'location', 'northeast')

    end
    title('Dry Season')
end

x3 = 0.27;
x4 = 0.83;
sb1.Position = [0.06    0.1207    x3    x4];
sb2.Position = [0.39    0.1207    x3    x4];
sb3.Position = [0.72    0.1207    x3    x4];

exportgraphics(f1, '../plots/amount_seasons.eps')

% subplot(1,2,2)
% ameff1 = ameff1./numb1;
% bar(bin1+db1/2,ameff1)
% grid on
% xlim([0 1])
% ylim([-4 0])
% legend('Waikiki', 'HIG', 'Lyon', 'Maunawili', 'Kailua', 'location', 'southwest')
% ylabel(['\delta^{18}O (',char(8240),')'])
% 
% xlabel('Rain rate percentile')
% % text(1, 7, 'b)', 'fontsize', 19)
% text(-9, 0, 'b)', 'fontsize', 19)





%% Seasonal effects
 

d1_avg = zeros(Nlocs,3);
dw_avg = zeros(Nlocs,3);
dd_avg = zeros(Nlocs,3);

d1_std = zeros(Nlocs,3);
dw_std = zeros(Nlocs,3);
dd_std = zeros(Nlocs,3);

d1_num = zeros(Nlocs,1);
dw_num = zeros(Nlocs,1);
dd_num = zeros(Nlocs,1);


for ns = 1:Nlocs
    
    
    mn1   = mn_all{ns};
    vec_d = mn1 > 4 & mn1 <= 10;
    vec_w = mn1 <= 4 | mn1 > 10;
    
    
    dD1   = dD_all{ns};
    d18o1 = d18O_all{ns};
    dexc1 = dexc_all{ns};
    rmm1  = rmm_all{ns};
    
    dDw   = dD1(vec_w);
    d18ow = d18o1(vec_w);
    dexcw = dexc1(vec_w);
    rmmw  = rmm1(vec_w);
    
    dDd   = dD1(vec_d);
    d18od = d18o1(vec_d);
    dexcd = dexc1(vec_d);
    rmmd  = rmm1(vec_d);
    
    
    % Total average
    vec   = ~isnan(dD1.*d18o1);
    d1_num(ns)   = nnz(vec);
    
    
    tmp1  = sum(dD1(vec).*rmm1(vec))./sum(rmm1(vec));
    d1_avg(ns,1) = tmp1;
    
    tmp2  = sum((dD1(vec)-tmp1).^2.*rmm1(vec))./sum(rmm1(vec));
    d1_std(ns,1) = sqrt(tmp2);
    
    tmp1  = sum(d18o1(vec).*rmm1(vec))./sum(rmm1(vec));
    d1_avg(ns,2) = tmp1; 
    
    tmp2  = sum((d18o1(vec)-tmp1).^2.*rmm1(vec))./sum(rmm1(vec));
    d1_std(ns,2) = sqrt(tmp2);
    
    tmp1  = sum(dexc1(vec).*rmm1(vec))./sum(rmm1(vec));
    d1_avg(ns,3) = tmp1; 
    
    tmp2  = sum((dexc1(vec)-tmp1).^2.*rmm1(vec))./sum(rmm1(vec));
    d1_std(ns,3) = sqrt(tmp2);
    
    
    
    % Dry season
    vec   = ~isnan(dDd.*d18od);
    dd_num(ns)   = nnz(vec);

    tmp1  = sum(dDd(vec).*rmmd(vec))./sum(rmmd(vec));
    dd_avg(ns,1) = tmp1;
    
    tmp2  = sum((dDd(vec)-tmp1).^2.*rmmd(vec))./sum(rmmd(vec));
    dd_std(ns,1) = sqrt(tmp2);
    
    tmp1  = sum(d18od(vec).*rmmd(vec))./sum(rmmd(vec));
    dd_avg(ns,2) = tmp1; 
    
    tmp2  = sum((d18od(vec)-tmp1).^2.*rmmd(vec))./sum(rmmd(vec));
    dd_std(ns,2) = sqrt(tmp2);

    tmp1  = sum(dexcd(vec).*rmmd(vec))./sum(rmmd(vec));
    dd_avg(ns,3) = tmp1; 
    
    tmp2  = sum((dexcd(vec)-tmp1).^2.*rmmd(vec))./sum(rmmd(vec));
    dd_std(ns,3) = sqrt(tmp2);
    
    
    
    % Wet season
    vec   = ~isnan(dDw.*d18ow);
    dw_num(ns)   = nnz(vec);
    
    tmp1  = sum(dDw(vec).*rmmw(vec))./sum(rmmw(vec));
    dw_avg(ns,1) = tmp1;

    tmp2  = sum((dDw(vec)-tmp1).^2.*rmmw(vec))./sum(rmmw(vec));
    dw_std(ns,1) = sqrt(tmp2);
    
    tmp1  = sum(d18ow(vec).*rmmw(vec))./sum(rmmw(vec));
    dw_avg(ns,2) = tmp1; 

    tmp2  = sum((d18ow(vec)-tmp1).^2.*rmmw(vec))./sum(rmmw(vec));
    dw_std(ns,2) = sqrt(tmp2);
    
    tmp1  = sum(dexcw(vec).*rmmw(vec))./sum(rmmw(vec));
    dw_avg(ns,3) = tmp1; 
   
    tmp2  = sum((dexcw(vec)-tmp1).^2.*rmmw(vec))./sum(rmmw(vec));
    dw_std(ns,3) = sqrt(tmp2);
    
end

%% LMWL by year 

v_18o0 = zeros(1,0);
v_D0   = zeros(1,0);

v_18o1 = zeros(1,0);
v_D1   = zeros(1,0);

v_18o2 = zeros(1,0);
v_D2   = zeros(1,0);

day0  = min(deploy);
day1  = day0+365;


close all
f = figure('Position', [1 1 600 600]);
hold on
for ns = Nlocs:-1:1 %1:Nlocs
    
    
    day   = day_all{ns};
    
    id1   = find(day >= day0 & day < day1);
    id2   = find(day >= day1);
    
    
    dD    = dD_all{ns};
    d18o  = d18O_all{ns};
    
    dD1   = dD(id1);
    d18o1 = d18o(id1);

    dD2   = dD(id2);
    d18o2 = d18o(id2);
    
    if(ns < 5)
        v_D0   = cat(2,v_D0,dD);
        v_D1   = cat(2,v_D1,dD1);
        v_D2   = cat(2,v_D2,dD2);
        v_18o0 = cat(2,v_18o0,d18o);
        v_18o1 = cat(2,v_18o1,d18o1);
        v_18o2 = cat(2,v_18o2,d18o2);
    else
        v_D0   = cat(2,v_D0,dD');
        v_D1   = cat(2,v_D1,dD1');
        v_D2   = cat(2,v_D2,dD2');
        v_18o0 = cat(2,v_18o0,d18o');
        v_18o1 = cat(2,v_18o1,d18o1');
        v_18o2 = cat(2,v_18o2,d18o2');
    end
    
    
    s(ns) = scatter(d18o1, dD1,  150, 'x', 'linewidth', 2);
%     s(ns).MarkerFaceColor = cmp(ns,:);
    s(ns).MarkerEdgeColor = cmp(ns,:);
    s(ns).MarkerFaceAlpha = 1;
    
    s2 = scatter(d18o2, dD2,  150, 'x', 'linewidth', 2);
%     s2.MarkerFaceColor = cmp(ns,:);
    s2.MarkerEdgeColor = cmp(ns,:);
    s2.MarkerFaceAlpha = 1;
    
end

grid on
box on
xlabel(['\delta^{18}O (',char(8240),')'])
ylabel(['\delta^{2}H (',char(8240),')'])


mdl_yr0 = fitlm(v_18o0,v_D0);
mdl_yr1 = fitlm(v_18o1,v_D1);
mdl_yr2 = fitlm(v_18o2,v_D2);

beta  = mdl_yr0.Coefficients.Estimate(1);
alpha = mdl_yr0.Coefficients.Estimate(2);

x0    = -10;
x1    = 10;

y0    = alpha*x0+beta;
y1    = alpha*x1+beta;

line([x0 x1],[y0 y1],  'color', 'k')

xlim([-10 6])
ylim([-80 40])

set(gca, 'xtick', -10:2:6)

legend([s(1) s(2) s(3) s(4) s(5)], 'Waikiki', 'HIG', 'Lyon', 'Maunawili', 'Kailua', 'location', 'southeast')

rectangle('Position', [-9.2 18 7.4 15], 'FaceColor', 'w')

text(-9.0, 30,'\delta^{2}H = 7.36*\delta^{18}O + 10.32', 'fontsize', 19)
text(-9.0, 22, 'R^{2} = 0.90', 'fontsize', 19)

exportgraphics(f, '../plots/lmwl.eps')





%% LMWL by season 

mdl_w = cell(1,Nlocs);
mdl_d = cell(1,Nlocs);

v_18ow = zeros(1,0);
v_Dw   = zeros(1,0);
v_18od = zeros(1,0);
v_Dd   = zeros(1,0);


close all
f = figure('Position', [1 1 600 600]);
hold on
for ns = 1:Nlocs
    
    
    mn1   = mn_all{ns};
    vec_d = mn1 > 4 & mn1 <= 10;
    vec_w = mn1 <= 4 | mn1 > 10;
    
    
    dD1   = dD_all{ns};
    d18o1 = d18O_all{ns};
    dexc1 = dexc_all{ns};
    rmm1  = rmm_all{ns};
    
    dDw   = dD1(vec_w);
    d18ow = d18o1(vec_w);
    dexcw = dexc1(vec_w);
    rmmw  = rmm1(vec_w);
    
    dDd   = dD1(vec_d);
    d18od = d18o1(vec_d);
    dexcd = dexc1(vec_d);
    rmmd  = rmm1(vec_d);

    if(ns < 5)
        v_Dw   = cat(2,v_Dw,dDw);
        v_Dd   = cat(2,v_Dd,dDd);
        v_18ow = cat(2,v_18ow,d18ow);
        v_18od = cat(2,v_18od,d18od);
    else
        v_Dw   = cat(2,v_Dw,dDw');
        v_Dd   = cat(2,v_Dd,dDd');
        v_18ow = cat(2,v_18ow,d18ow');
        v_18od = cat(2,v_18od,d18od');
    end
    
    s1 = scatter(dDw, d18ow, 150, 'c');
    s1.MarkerFaceColor = cmp(ns,:);
    s1.MarkerEdgeColor = 'k';

    
    s2 = scatter(dDd, d18od, 150, 'd');
    s2.MarkerFaceColor = cmp(ns,:);
    s2.MarkerEdgeColor = 'k';
    
    
    mdl_w{ns} = fitlm(d18ow,dDw);
    mdl_d{ns} = fitlm(d18od,dDd);
    
    
end

grid on
box on
ylabel(['\delta^{18}O (',char(8240),')'])
xlabel(['\delta^{2}H (',char(8240),')'])


mdl_tot_w = fitlm(v_18ow,v_Dw);
mdl_tot_d = fitlm(v_18od,v_Dd);

exportgraphics(f, '../plots/lmwl_season.eps')


%% LMWL by location 


v_18o0 = zeros(1,0);
v_D0   = zeros(1,0);

v_18o1 = zeros(1,0);
v_D1   = zeros(1,0);

v_18o2 = zeros(1,0);
v_D2   = zeros(1,0);


close all
for ns = 1:Nlocs
        
    dD    = dD_all{ns};
    d18o  = d18O_all{ns};
    
    if(ns <= 3)
        v_D0   = cat(2, v_D0,   dD);
        v_18o0 = cat(2, v_18o0, d18o);
    elseif(ns == 4)
        v_D1   = cat(2, v_D1,   dD);
        v_18o1 = cat(2, v_18o1, d18o);
    elseif(ns == 5)
        v_D1   = cat(2, v_D1,   dD');
        v_18o1 = cat(2, v_18o1, d18o');
        
    end
    
    
end


mdl_lee  = fitlm(v_18o0,v_D0);
mdl_wind  = fitlm(v_18o1,v_D1);


%% Misc


%%





day     = mat(:,3)+d0-d1-1;

mn      = month(datetime(datestr(day)));


vec_l   = ~isnan(dD);

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


%% Compare rainfall abundances by season


day0 = '01-Jan-2019';

day1 = '01-Apr-2019';
day2 = '01-Oct-2019';
day3 = '01-Apr-2020';
day4 = '01-Oct-2020';
day5 = '01-Apr-2021';
day6 = '01-Oct-2022';

d0   = datenum(day0);
d1   = datenum(day1)-d0+1;
d2   = datenum(day2)-d0+1;
d3   = datenum(day3)-d0+1;
d4   = datenum(day4)-d0+1;
d5   = datenum(day5)-d0+1;
d6   = datenum(day6)-d0+1;


db1       = 1;
bin1      = -100:db1:100;
Nb1       = length(bin1);

db2       = 0.3;
bin2      = -100:db2:100;
Nb2       = length(bin2);


hc1d      = zeros(Nb1-1,Nlocs);
hc2d      = zeros(Nb2-1,Nlocs);
hc3d      = zeros(Nb1-1,Nlocs);

hc1w      = zeros(Nb1-1,Nlocs);
hc2w      = zeros(Nb2-1,Nlocs);
hc3w      = zeros(Nb1-1,Nlocs);


for ns = 1:Nlocs

    vec0 = cell2mat(day_all(ns));
    vec1 = cell2mat(dD_all(ns));
    vec2 = cell2mat(d18O_all(ns));
    vec3 = cell2mat(dexc_all(ns));

    ids  = find(~isnan(vec1) & ((vec0 >= d1 & vec0 < d2) |  ...
        (vec0 >= d3 & vec0 < d4) | ...
        (vec0 >= d5 & vec0 < d6)) );
    
    idw  = find(~isnan(vec1) & ~((vec0 >= d1 & vec0 < d2) |  ...
        (vec0 >= d3 & vec0 < d4) | ...
        (vec0 >= d5 & vec0 < d6)) );
    
    hc1d(:,ns) = histcounts(vec1(ids), bin1, 'normalization', 'cdf');
    hc2d(:,ns) = histcounts(vec2(ids), bin2, 'normalization', 'cdf');
    hc3d(:,ns) = histcounts(vec3(ids), bin1, 'normalization', 'cdf');
    
    hc1w(:,ns) = histcounts(vec1(idw), bin1, 'normalization', 'cdf');
    hc2w(:,ns) = histcounts(vec2(idw), bin2, 'normalization', 'cdf');
    hc3w(:,ns) = histcounts(vec3(idw), bin1, 'normalization', 'cdf');
    
end


%% 

close all
f = figure('Position', [1 1 1200 700]);
subplot(2,2,1)
plot(bin1(1:end-1), hc1w, '-o')
grid on
ylim([0 1])
xlim([-15 10])
ylabel('CDF')
xlabel(['\delta^{2}H (',char(8240),')'])
text(-14.6, 0.94, 'a)', 'fontsize', 19)
title('Wet Season')

subplot(2,2,2)
plot(bin1(1:end-1), hc1d, '-o')
grid on
ylim([0 1])
xlim([-15 10])
ylabel('CDF')
xlabel(['\delta^{2}H (',char(8240),')'])
text(-14.6, 0.94, 'c)', 'fontsize', 19)
title('Dry Season')

subplot(2,2,3)
plot(bin2(1:end-1), hc2w, '-d')
grid on
ylim([0 1])
xlim([-5 0])
ylabel('CDF')
xlabel(['\delta^{18}O (',char(8240),')'])
text(-4.87, 0.94, 'b)', 'fontsize', 19)
legend('Waikiki', 'HIG', 'Lyon', 'Maunawili', 'Kailua', 'location', 'southeast')

subplot(2,2,4)
plot(bin2(1:end-1), hc2d, '-d')
grid on
ylim([0 1])
xlim([-5 0])
ylabel('CDF')
xlabel(['\delta^{18}O (',char(8240),')'])
text(-4.87, 0.94, 'd)', 'fontsize', 19)
exportgraphics(f, '../plots/seasonal_abundances_cdf.eps')


%%

close all
f = figure('Position', [1 1 1200 400]);
subplot(1,2,1)
plot(bin1(1:end-1), hc3w, '-d')
grid on
ylim([0 1])
xlim([0 20])
set(gca, 'xtick', 0:4:20)
ylabel('CDF')
xlabel(['d-excess (',char(8240),')'])
text(0.3, 0.97, 'a)', 'fontsize', 19)
title('Wet Season')

subplot(1,2,2)
plot(bin1(1:end-1), hc3d, '-d')
grid on
ylim([0 1])
xlim([0 20])
set(gca, 'xtick', 0:4:20)
ylabel('CDF')
xlabel(['d-excess (',char(8240),')'])
text(0.3, 0.97, 'b)', 'fontsize', 19)
title('Dry Sesason')
legend('Waikiki', 'HIG', 'Lyon', 'Maunawili', 'Kailua', 'location', 'southeast')
exportgraphics(f, '../plots/seasonal_dexc_cdf.eps')

