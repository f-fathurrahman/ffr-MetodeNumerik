%IT_SOLVE_STEP_NAVIER iterative solution of N-S in step domain
%   IFISS scriptfile: HCE;  DJS; 2 October 2013
% Copyright (c) 2013 D.J. Silvester, H.C. Elman, A. Ramage
clear variables
fprintf('FAST solution of flow problem in step domain ...\n')
viscosity=default('viscosity parameter (default 1/50)',1/50);
nlmethod=default('Picard/Newton/hybrid linearization 1/2/3 (default hybrid)',3);
nlmethod=nlmethod-1;
if nlmethod==0,
   maxit_p=default('number of Picard iterations (default 9)',9);
   maxit_n=0;
elseif nlmethod==1,
   maxit_p=0;
   maxit_n=default('number of Newton iterations (default 6)',6);
else
   maxit_p=default('number of Picard iterations (default 2)',2);
   maxit_n=default('number of Newton iterations (default 4)',4);
end
tol_nl=default('nonlinear tolerance (default 1.e-5)',1.e-5);
%
% set indicator for strategy for stopping criterion for Picard/Newton steps
exact = default('inexact or exact linear iteration 1/2 (default inexact)',1); 
%
% set parameters for iterative solution algorithm
ns_params.itmeth = 1;         % GMRES/Bicgstab(2)  1/2
ns_params.tol    = 1.d-6;     % For stokes only, set below for NS
ns_params.maxit  = 200;
ns_params.precon = 4;         % no / BFBT / Fp / xBFBt / PCD* / LSC* 0/1/2/3/4/5
ns_params.precon_format = 1;  % ideal / amg  preconditioning? 1/2
%
%% initialize for nonlinear iteration: compute Stokes solution
%% load assembled matrices
gohome
cd datafiles
load step_stokes_nobc.mat
%geometry_data = struct('domain',domain, 'x',x, 'y',y,'xy',xy,'xyp',xyp, ...
%                       'mv',mv,'mp',mp, 'outbnd',outbnd, 'qmethod',qmethod);
%
fprintf('stokes system ...\n') 
%% boundary conditions
[Ast,Bst,fst,gst] = flowbc(A,B,f,g,xy,bound);
nlres0_norm = norm([fst;gst]);
%
nv=length(fst)/2; np=length(gst); 
if qmethod>1,
   beta=0;
%  xst=[Ast,Bst';Bst,sparse(np,np)]\[fst;gst];
   xst=it_nstokes(Ast,Bst,sparse(np,np),beta,Q,G, [fst;gst], ...
                  zeros(length(Ast),1),1,ns_params, domain); 
elseif qmethod==1
   beta=1/4;     % default parameter
%  xst=[Ast,Bst';Bst,-beta*C]\[fst;gst];
   xst=it_nstokes(Ast,Bst,C,beta,Q,G, [fst;gst], ...
                  zeros(length(Ast),1),1,ns_params, domain);
elseif qmethod==0
   fprintf('computing pressure stabilized solution...\n')
   beta=1;
%  xst=[Ast,Bst';Bst,-C]\[fst;gst]; beta=1;
   xst=it_nstokes(Ast,Bst,C,beta,Q,G, [fst;gst], ...
                  zeros(length(Ast),1),1,ns_params, domain); 
end
        xstit=xst;
% compute residual of Stokes solution
if qmethod>1
   N = navier_q2(xy,mv,xst);
elseif qmethod<=1,
   nubeta=beta/viscosity;
fprintf('\n nubeta is set to %e \n',nubeta)
   N = navier_q1(xy,ev,xst);
end
Anst = viscosity*A + [N, sparse(nv,nv); sparse(nv,nv), N];
[Anst,Bst,fst,gst] = flowbc(Anst,B,f,g,xy,bound);
if     qmethod>1,  nlres = [Anst,Bst';Bst,sparse(np,np)]*xst-[fst;gst];
elseif qmethod<=1, nlres = [Anst,Bst';Bst,-nubeta*C]*xst-[fst;gst];
end
nlres_norm  = norm(nlres);
%%% plot solution
spc=default('uniform/nonlinear streamlines 1/2 (default uniform)',1);
contourn = default('number of contour lines (default 50)',50);
flowplot09(qmethod,xst,By,Bx,A,xy,xyp,x,y,bound,bndxy,bnde,obs,contourn,spc,33)
%flowplotl(qmethod,xst,By,Bx,A,xy,xyp,x,y,bound,33);
pause(1) 
fprintf('\n\ninitial nonlinear residual is %e ',nlres0_norm)
fprintf('\nStokes solution residual is %e\n', nlres_norm)
flowsol = xst;
%
%
pde=4;
it_p = 0;
tau_safety = 1.d-2;
%
% nonlinear iteration 
%% Picard startup step
while nlres_norm>nlres0_norm*tol_nl && it_p<maxit_p,
   it_p = it_p+1;
   fprintf('\nPicard iteration number %g \n',it_p),
   % stopping tolerance
   if exact~=1, ns_params.tol = tol_nl/10;
   else         ns_params.tol = tau_safety*nlres_norm;
   end
%
% compute Picard correction and update solution
  %if     qmethod>1,  dxns = -[Anst,Bst';Bst,sparse(np,np)]\nlres;
   if qmethod>1, 
       dxns = -it_nstokes(Anst,Bst,sparse(np,np),beta,Q,G, nlres, ...
                          flowsol,viscosity,ns_params, domain);
  %elseif qmethod<=1, dxns = -[Anst,Bst';Bst,-nubeta*C]\nlres;
   elseif qmethod<=1, 
      dxns = -it_nstokes(Anst,Bst,C,nubeta,Q,G, nlres, ...
                         flowsol,viscosity,ns_params, domain);
   end
   xns = flowsol + dxns;
% compute residual of new solution
   if     qmethod>1,  N = navier_q2(xy,mv,xns);
   elseif qmethod<=1, N = navier_q1(xy,ev,xns);
   end
   Anst = viscosity*A + [N, sparse(nv,nv); sparse(nv,nv), N];
   [Anst,Bst,fst,gst] = flowbc(Anst,B,f,g,xy,bound);
   if     qmethod>1,  nlres = [Anst,Bst';Bst,sparse(np,np)]*xns-[fst;gst];
   elseif qmethod<=1, nlres = [Anst,Bst';Bst,-nubeta*C]*xns-[fst;gst];
   end
   nlres_norm = norm(nlres);
   nnv=length(fst); soldiff=norm(xns(1:nnv)-flowsol(1:nnv));
   fprintf('nonlinear residual is %e',nlres_norm)
   fprintf('\n   velocity change is %e\n',soldiff)
% plot solution
   flowplot09(qmethod,xns,By,Bx,A,xy,xyp,x,y,bound,bndxy,bnde,obs,contourn,spc,66);drawnow
%  flowplotz(qmethod,xns,By,Bx,A,xy,xyp,x,y,bound,66); drawnow
   pause(1)
   flowsol = xns;
% end of Picard iteration loop
end
%%
%
it_nl = it_p;
it_n = 0;
%% Newton iteration loop
while (nlres_norm > nlres0_norm*tol_nl) && (it_nl < maxit_p + maxit_n),
   it_n = it_n+1;
   it_nl = it_nl+1;
   fprintf('\nNewton iteration number %g \n',it_n),
% stopping tolerance
   if exact~=1, ns_params.tol = tol_nl/10;
   else         ns_params.tol = tau_safety*nlres_norm;
   end
% compute Jacobian of current solution
   if     qmethod>1,  [Nxx,Nxy,Nyx,Nyy] = newton_q2(xy,mv,flowsol);
   elseif qmethod<=1, [Nxx,Nxy,Nyx,Nyy] = newton_q1(xy,ev,flowsol);
   end
   J = viscosity*A + [N + Nxx, Nxy; Nyx, N + Nyy];
   Jnst = newtonbc(J,xy,bound); 
% compute Newton correction and update solution
  %if     qmethod>1,  dxns = -[Jnst,Bst';Bst,sparse(np,np)]\nlres;
   if qmethod>1,  
      dxns = -it_nstokes(Jnst,Bst,sparse(np,np),beta,Q,G, nlres,...
                         flowsol,viscosity,ns_params, domain);
  %elseif qmethod<=1, dxns = -[Jnst,Bst';Bst,-nubeta*C]\nlres;
   elseif qmethod<=1, 
      dxns = -it_nstokes(Jnst,Bst,C,nubeta,Q,G, nlres, ...
                         flowsol,viscosity,ns_params, domain);
   end
   xns = flowsol + dxns;
% compute residual of new solution
   if     qmethod>1,  N = navier_q2(xy,mv,xns);
   elseif qmethod<=1, N = navier_q1(xy,ev,xns);
   end
   Anst = viscosity*A + [N, sparse(nv,nv); sparse(nv,nv), N];
   [Anst,Bst,fst,gst] = flowbc(Anst,B,f,g,xy,bound);
   if     qmethod>1,  nlres = [Anst,Bst';Bst,sparse(np,np)]*xns-[fst;gst];
   elseif qmethod<=1, nlres = [Anst,Bst';Bst,-nubeta*C]*xns-[fst;gst];
   end
   nlres_norm = norm(nlres);
   nnv=length(fst); soldiff=norm(xns(1:nnv)-flowsol(1:nnv));
   fprintf('nonlinear residual is %e',nlres_norm)
   fprintf('\n   velocity change is %e\n',soldiff)
% plot solution
   flowplot09(qmethod,xns,By,Bx,A,xy,xyp,x,y,bound,bndxy,bnde,obs,contourn,spc,66);drawnow
%   flowplotz(qmethod,xns,By,Bx,A,xy,xyp,x,y,bound,66); drawnow
   pause(1)
   flowsol = xns;
% end of Newton iteration loop 
end
%%
if nlres_norm <= nlres0_norm * tol_nl, 
   fprintf('\nfinished, nonlinear convergence test satisfied\n\n');
else
   fprintf('\nfinished, stopped on iteration counts\n\n');
end
%
xitns=flowsol;
%% estimate errors
if qmethod==1
   [jmpx,jmpy,els] = stressjmps_q1p0(viscosity,flowsol,xy,ev,ebound);
   [error_x,error_y,fex,fey,ae] = navierpost_q1p0_p(viscosity,flowsol,jmpx,jmpy,els,xy,ev);
   [error_x,error_y] = navierpost_q1p0_bc(viscosity,ae,fex,fey,...
                                        error_x,error_y,xy,ev,ebound);
   error_div = q1div(xy,ev,flowsol);
   errorestit=sqrt(sum(error_x.^2 + error_y.^2 + error_div.^2));
   fprintf('estimated overall error is %10.6e \n',errorestit)
   ee_errorit=sqrt((error_x.^2 + error_y.^2 + error_div.^2));
%% plot element errors
   %eplotl(ee_error,ev,xy,x,y,67);
%% plot macroelement errors
   mplotl(ee_errorit,ev,xy,x,y,67), title('Estimated error')
   pause(5), figure(66)
elseif qmethod==0
   error_div = q1div(xy,ev,flowsol);
elseif qmethod>1, navierpost, 	
end
