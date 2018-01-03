x = 0;
y = [7.15014e6 0 0 0.937045e-3];
xStop = 1200;
h = 50;
freq = 2;
[xSol,ySol] = runKut4(@fex7_6,x,y,xStop,h);
printSol(xSol,ySol,freq)