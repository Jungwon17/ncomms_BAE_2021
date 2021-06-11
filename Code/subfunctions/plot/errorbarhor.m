function h = errorbarhor(data, yP, barWidth)
nData = sum(~isnan(data));
mR = nanmean(data);
sR = nanstd(data) / sqrt(nData);

x = [mR-sR, mR, mR+sR, NaN, mR-sR, mR-sR, NaN, mR+sR, mR+sR];
y = [yP, yP, yP, NaN, yP-barWidth, yP+barWidth, NaN, yP-barWidth, yP+barWidth];

h(1) = plot(x, y, 'LineWidth', 0.35, 'Color', [0 0 0]);
h(2) = plot(x(2), y(2), 'Marker', 'o', 'MarkerSize', 1, 'Color', [0 0 0]);