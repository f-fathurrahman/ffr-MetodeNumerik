% Example 3.6 (Interpolation)
xData = [0.1; 0.2; 0.5; 0.6; 0.8; 1.2; 1.5];
yData = [-1.5342; -1.0811; -0.4445; -0.3085; -0.0868; 0.2281; 0.3824];

x = 0.1:0.05:1.5;

n = length(x);
y = zeros(n,2);

for i = 1:n
  y(i,1) = rational_interp(xData,yData,x(i));
  y(i,2) = neville(xData,yData,x(i));
end

plot(x,y(:,1),'k-'); hold on
plot(x,y(:,2),'k:'); hold on

plot(xData, yData, 'ko')

grid on
xlabel('x');
ylabel('y')