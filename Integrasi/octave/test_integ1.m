% Example 6.4 (Recursive trapezoidal rule)

func = @(x) (sqrt(x)*cos(x));

I2h = 0;
for k = 1:20
  Ih = integ_trapezoid_recursive(func,0,pi,I2h,k);
  diffI = abs(Ih - I2h);
  fprintf('k = %5d %18.10f %18.10e\n', k, Ih, diffI);
  if( k > 1 && diffI < 1.0e-6 )
    Integral = Ih;
    No_of_func_evaluations = 2^(k-1) + 1;
    return
  end
  I2h = Ih;
end
error('Too many iterations')


