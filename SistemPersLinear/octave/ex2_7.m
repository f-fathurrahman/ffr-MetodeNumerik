% Example 2.7 (Doolittle’s decomposition)
A = [3 -1 4; -2 0 5; 7 2 -2];
B = [6 -4; 3 2; 7 -5];
A = LUdec(A);
det1 = prod(diag(A))
for i = 1:size(B,2)
  X(:,i) = LUsol(A,B(:,i));
end
X