% Needed for gaussPiv
addpath('../../SistemPersLinear/octave')

% Example 3.12 (Polynomial curve fitting)
xData = [-0.04,0.93,1.95,2.90,3.83,5.0, 5.98,7.05,8.21,9.08,10.09]';
yData = [-8.66,-6.44,-4.36,-3.27,-0.88,0.87,3.31,4.63,6.19,7.4,8.85]';

format short e
while 1
  k = input('degree of polynomial = ');
  if isempty(k)
    fprintf('Done\n\n')
    % Loop is terminated
    % by pressing ''return''
    break
  end
  coeff = polynFit(xData,yData,k+1)
  sigma = stdDev(coeff,xData,yData)
  fprintf('\n')
end
