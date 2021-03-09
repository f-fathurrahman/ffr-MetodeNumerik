x = 0:0.01:1;
y = 0:0.01:2;

F = zeros(length(x), length(y));
for i = 1:length(x)
  for j = 1:length(y)
    F(i,j) = x(i)^2 + y(j)^2;
  end
end

p.x = x;
p.y = y;
p.dim = 2;
p.measure = 2;
data.points = p;
data.values = F;

disp('Integral: ')
v1 = riemann.trpzd(data)
disp('Mean value:')
v2 = riemann.mean_value(data)
