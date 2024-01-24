% Example 7.11 (Bulirsch-Stoer integration)
[xSol,ySol] = bulStoer(@fex7_11, 0, [0 0], 10, 0.5);
printSol(xSol,ySol,1)

plot(xSol,ySol(:,2),'k-o')
grid on
xlabel('Time (s)')
ylabel('Current (A)')