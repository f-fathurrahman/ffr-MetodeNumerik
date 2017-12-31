% Example 10.7 (Minimization with downhill simplex)
mu = 10000;
S = @(x) x(1) + 2*x(2)*sec(x(3)) + mu*((x(1) + x(2)*tan(x(3)))*x(2) - 8)^2;
xStart = [4;2;0];
[x,fMin,nCyc] = downhill(S,xStart);
perim = x(1) + 2*x(2)*sec(x(3));
area = (x(1) + x(2)*tan(x(3)))*x(2);
fprintf('b = %8.4f\n',x(1))
fprintf('h = %8.4f\n',x(2))
fprintf('theta_in_deg = %8.4f\n',x(3)*180/pi)
fprintf('perimeter = %8.4f\n',perim)
fprintf('area = %8.4f\n',area)
fprintf('number of cycles = %2.0f\n',nCyc)