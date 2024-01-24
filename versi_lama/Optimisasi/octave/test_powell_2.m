% Example 10.4 (Powell's method of minimization)
mu = 1.0;

func = @(x) (x(1) - 5.0)^2 + (x(2) - 8.0)^2 + mu*(x(1)*x(2) - 5.0)^2;

xStart = [1.0;5.0];

[x,f,nCyc] = powell(func,xStart);

fprintf('Intersection point = %8.5f %8.5f\n',x(1),x(2))

xy = x(1)*x(2);
fprintf('Constraint x*y = %8.5f\n',xy)

dist = sqrt((x(1) - 5.0)^2 + (x(2) - 8.0)^2);
fprintf('Distance = %8.5f\n',dist)

fprintf('Number of cycles = %2.0f\n',nCyc)