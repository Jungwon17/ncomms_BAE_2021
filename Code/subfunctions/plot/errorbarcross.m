function h = errorbarcross(dataX, dataY)
nDataX = sum(~isnan(dataX));
mX = nanmean(dataX);
sX = nanstd(dataX) / sqrt(nDataX);

nDataY = sum(~isnan(dataY));
mY = nanmean(dataY);
sY = nanstd(dataY) / sqrt(nDataY);

x = [mX-sX, mX+sX, NaN, mX, mX];
y = [mY, mY, NaN, mY-sY, mY+sY];

% h(1) = plot(mX, mY, 'Marker', 'o', 'MarkerSize', 1, 'Color', [0 0 0]);
h = plot(x, y, 'LineWidth', 0.35, 'Color', [0 0 0]);
