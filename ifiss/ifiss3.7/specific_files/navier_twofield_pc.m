%NAVIER_TWOFIELD_PC sets up two field pressure Navier-Stokes test cases
%   IFISS scriptfile: DJS; 5 August 2023.
% Copyright (c) 2005 D.J. Silvester, H.C. Elman, A. Ramage 
gohome
clear variables
fprintf('\nspecification of two-field Navier-Stokes problem.\n')
fprintf('\nchoose specific example (default is cavity)');
%fprintf('\n     1  Channel domain')
fprintf('\n     2  Flow over a backward facing step')
fprintf('\n     3  Lid driven cavity\n')
sn = default('',3);
if sn==2
system('copy .\stokes_flow\test_problems\backwardstep_flow.m .\stokes_flow\specific_flow.m');
system('copy .\stokes_flow\test_problems\backwardstep_bc.m .\stokes_flow\stream_bc.m');
   stepX_stokes, solveX_step_navier
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
   squareX_stokes, solveX_navier
else
   error('reference problem datafile not found!')
end
