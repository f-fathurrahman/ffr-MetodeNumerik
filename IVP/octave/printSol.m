function printSol(xSol,ySol,freq)
% Prints xSol and ySoln arrays in tabular format.
% USAGE: printSol(xSol,ySol,freq)
% freq = printout frequency (prints every freq-th
%        line of xSol and ySol).

[m,n] = size(ySol);

if freq == 0
  freq = m;
end

head = '     x';

for i = 1:n
  head = strcat(head,'      y',num2str(i));
end

fprintf(head)
fprintf('\n')

for i = 1:freq:m
  fprintf('%14.4e',xSol(i),ySol(i,:));
  fprintf('\n')
end

if i ~= m
  fprintf('%14.4e',xSol(m),ySol(m,:));
end
