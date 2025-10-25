%CD_TESTPROBLEM_UNIX sets up Examples 6.1 to 6.4
%   IFISS scriptfile: DJS; 9 May 2012;  19 September 2019
% Copyright (c) 2005 D.J. Silvester, H.C. Elman, A. Ramage 
gohome
clear variables
fprintf('\nspecification of reference convection-diffusion problem.\n')
fprintf('\nchoose specific example');
fprintf('\n     1  Constant vertical wind')
fprintf('\n     2  Vertical wind, characteristic layers')
fprintf('\n     3  Constant wind @ 30 degree angle')
fprintf('\n     4  Recirculating wind\n');
sn = default('',1);

if sn==3 
system('/bin/cp ./convection/test_problems/constant_wind.m ./convection/specific_wind.m');
system('/bin/cp ./convection/test_problems/test_bc.m ./diffusion/specific_bc.m');
   square_cd, solve_cd
elseif sn==4
system('/bin/cp ./convection/test_problems/circular_wind.m ./convection/specific_wind.m');
hotwall_model=default('regularise hot wall corner singularities? 1/0 (no)',0);
if hotwall_model==1
system('/bin/cp ./convection/test_problems/hotwallx_bc.m ./diffusion/specific_bc.m');
else
system('/bin/cp ./convection/test_problems/hotwall_bc.m ./diffusion/specific_bc.m');
end
   square_cd, solve_cd
elseif sn==1
system('/bin/cp ./convection/test_problems/vertical_wind.m ./convection/specific_wind.m');
system('/bin/cp ./convection/test_problems/solution1_bc.m ./diffusion/specific_bc.m');
   ref_cd, solve_cd
elseif sn==2
system('/bin/cp ./convection/test_problems/ref_wind.m ./convection/specific_wind.m');
system('/bin/cp ./convection/test_problems/ref_bc.m ./diffusion/specific_bc.m');
   ref_cd, solve_cd
else
   error('reference problem datafile not found!')
end
