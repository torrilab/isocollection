


v1 = [21.401619 21.37291 21.3330 21.2983 21.284012];
v2 = [-157.740208 -157.76576 -157.8025 -157.8166 -157.851432];
v3 = 300*ones(1,5);
v4 = [5 4 3 2 1];
% v4 = [cmp(1,:), cmp(2,:), cmp(3,:), cmp(4,:), cmp(5,:)];
set(0, 'defaultaxesfontsize', 19)
set(0, 'defaultlinelinewidth', 2)

close all

f = figure('Position', [1 1 900 700]);

g = geoscatter(v1,v2, v3, v4, 'fill');
g.MarkerEdgeColor = 'k';

geolimits([21.225 21.73],[-158.25 -157.65])
geobasemap topographic

cmp = colormap(lines(5));

set(gca, 'Fontsize', 19)

% legend('Kailua', 'Maunawili', 'Lyon', 'HIG', 'Lyon')

% exportgraphics(f, '../plots/oahu_map.eps')

