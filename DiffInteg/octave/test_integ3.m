% Function to be integrated, using change of variable
func = @(t) (2* t^2 * cos(t^2));

[Integral, numEval] = integ_romberg( func, 0, sqrt(pi) );

fprintf('Integral = %18.10f\n', Integral)
fprintf('numEval  = %d\n', numEval)


