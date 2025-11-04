%STOKES_TWOFIELD_PC sets up two field pressure Stokes test cases
%   IFISS scriptfile: DJS; 5 August 2023.
% Copyright (c) 2005 D.J. Silvester, H.C. Elman, A. Ramage 
gohome
clear variables
fprintf('\nspecification of two-field Stokes problem.\n')
fprintf('\nchoose specific example (default is cavity)');
%fprintf('\n     1  Channel domain')
fprintf('\n     2  Flow over a backward facing step')
fprintf('\n     3  Lid driven cavity')
fprintf('\n     4  Colliding flow\n');
sn = default('',3);
if sn==1, 
system('copy .\stokes_flow\test_problems\poiseuille_flow.m .\stokes_flow\specific_flow.m');
system('copy .\stokes_flow\test_problems\poiseuille_bc.m .\stokes_flow\stream_bc.m');
%   channel_stokes, xsolve_stokes
elseif sn==2 
system('copy .\stokes_flow\test_problems\backwardstep_flow.m .\stokes_flow\specific_flow.m');
system('copy .\stokes_flow\test_problems\backwardstep_bc.m .\stokes_flow\stream_bc.m');
   stepX_stokes, solveX_step_stokes
elseif sn==3
   lid_model=default('cavity type leaky/tight/regularised 1/2/3 (regularised)',3);
   if lid_model ==3    
system('copy .\stokes_flow\test_problems\regcavity_flow.m .\stokes_flow\specific_flow.m');
   elseif lid_model ==2    
system('copy .\stokes_flow\test_problems\tightcavity_flow.m .\stokes_flow\specific_flow.m');
   else
system('copy .\stokes_flow\test_problems\leakycavity_flow.m .\stokes_flow\specific_flow.m');
   end
system('copy .\stokes_flow\test_problems\zero_bc.m .\stokes_flow\stream_bc.m');
   squareX_stokes, solveX_stokes
elseif sn==4
system('copy .\stokes_flow\test_problems\collide_flow.m .\stokes_flow\specific_flow.m');
system('copy .\stokes_flow\test_problems\collide_bc.m .\stokes_flow\stream_bc.m');
   squareX_stokes, solveX_stokes
else
   error('reference problem datafile not found!')
end
