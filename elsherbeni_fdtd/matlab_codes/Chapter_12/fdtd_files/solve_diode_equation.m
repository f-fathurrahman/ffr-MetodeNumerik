function [x] = solve_diode_equation(A, B, C, x)
% Function used to solve the diode equation 
% which is in the form Ae^{Bx}+x+C=0
% using the Newton-Raphson method

tolerance = 1e-25;
max_iter = 50;
iter = 0;
f = A * exp(B*x) + x + C;
while ((iter < max_iter) && (abs(f) > tolerance))
    fp = A * B * exp(B*x) + 1;
    x = x - f/fp;
    f = A * exp(B*x) + x + C;
    iter = iter + 1;
end