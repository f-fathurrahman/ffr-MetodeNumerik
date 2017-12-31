% Example 3.4 (rational function interpolation)
xData = [0.15; 2.3; 3.15; 4.85; 6.25; 7.95];
yData = [4.79867; 4.49013; 4.22430; 3.47313; 2.66674; 1.51909];
for x = 0:0.5:8
  y = rational_interp(xData, yData, x);
  yExact = 4.8*cos(pi*x/20);
  fprintf('%10.5f %10.5f %10.5f    %10.5e\n', x, y, yExact, abs(y-yExact) );
end
