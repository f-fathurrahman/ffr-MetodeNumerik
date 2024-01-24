n = 20;
x = zeros(n,1);
b = zeros(n,1); b(n) = 1;
[x,numIter] = conjGrad(@fex2_18,x,b)