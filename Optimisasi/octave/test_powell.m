% Example 10.3 (Powell’s method of minimization)
func = @(x) 100.0*(x(2) - x(1)^2)^2 + (1.0 -x(1))^2;

[xMin,fMin,numCycles] = powell(func,[-1,1])
