%%ifissdemo  Boussinesq unsteady flow over a heated step 
global Ra Pr H L
close all
fprintf('\nBoussinesq flow over a heated step  ... \n'),
fprintf('running STABTR for 200 timesteps  ... \n'), pause(5)
batchmode('B-NS2'), load unsteadyrun.mat
figure(13), set(gcf,'Position',[900,450,350,350],'Visible','on');
fprintf('\nCHECK OUT the time step history and final solution \n'), pause(5)
figure(13), set(gcf,'Visible','off');
if isunix
fprintf('CHECK OUT the temperature evolution  \n'),
status=system('open TemperatureMovie.mp4'); pause
end
ppbouss_checkpointdata
save unsteadyrun.mat
fprintf('CHECK the iterative solver convergence ...\n'),
pause(5)
% iterative solvers
batchmode('snapshot_boussx1')
figure(19), set(gcf,'Position',[900,400,350,350]); drawnow, 
batchmode('snapshot_boussx2')
title('inexact preconditioning','FontSize',12),
legend( 'AMG-PCD*','AMG-LSC')
figure(19), set(gcf,'Position',[900,400,350,350]); drawnow, 
fprintf('End of Boussinesq flow demo. Voila!\n'),
