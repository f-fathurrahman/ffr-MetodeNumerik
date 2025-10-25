%SOLVEX_NAVIER solve two-field Navier-Stokes problem
%   IFISS scriptfile:  DJS; 4 August 2023
% Copyright (c) 2023 D.J. Silvester, J. Pestana
clear variables
fprintf('Enclosed flow problem ...\n')
viscosity=default('viscosity parameter (default 1/100)',1/100);
nlmethod=0;
   maxit_p=default('number of Picard iterations (default 5)',5);
   maxit_n=0;
tol_nl=default('nonlinear tolerance (default 10*eps)',10*eps);
%
%
%% initialize for nonlinear iteration: compute Stokes solution
%% load assembled matrices
gohome
cd datafiles
load squareX_stokes_nobc.mat
%
fprintf('stokes system ...\n')
B=[B1;B0]; g=[g1;g0];
np1=length(g1); np0=length(g0); nu=length(f)/2; np=np1+np0;
%-------------------- debug test
%B=B1; g=g1; np=np1;

%% boundary conditions
[Ast,Bst,fst,gst] = flowbc(A,B,f,g,xy,bound);
nlres0_norm = norm([fst;gst]);
%
nv=length(fst)/2;
%----- compute flow solution
np=length(gst); tic
%-------------------------------------------------------
% xst=[Ast,Bst';Bst,sparse(np,np)]\[fst;gst];
% stable solve
       tt1=[ones(np1,1);zeros(np0,1)]/np1;
       tt0=[zeros(np1,1);ones(np0,1)]/np0;
xx = [Ast,Bst',sparse(2*nu,2); ...
      Bst,sparse(np,np),tt1,tt0; ...
      sparse(1,2*nu),tt1',zeros(1,2); ...
      sparse(1,2*nu),tt0',zeros(1,2)]\[fst;gst;0;0];
      xst=xx(1:end-2); multiplier=xx(end-1:end); etoc=toc;
etoc=toc; fprintf('Stokes system solved in %8.3e seconds\n\n',etoc)
%
%----------- scale divergence contribution matrices
p0alpha = default('scaling parameter (viscosity)',viscosity);
fprintf('scaling approximation with %e ...\n\n',p0alpha)
B=[B1;p0alpha*B0]; g=[g1;p0alpha*g0];
   
% compute residual of Stokes solution
      N = navier_q2(xy,mv,xst);
%
Anst = viscosity*A + [N, sparse(nv,nv); sparse(nv,nv), N];
[Anst,Bst,fst,gst] = flowbc(Anst,B,f,g,xy,bound);
np=length(gst); nuv=length(fst);
nlres = [Anst,Bst';Bst,sparse(np,np)]*xst-[fst;gst];
nlres_norm  = norm(nlres);
fprintf('\n\ninitial nonlinear residual is %e ',nlres0_norm)
fprintf('\nStokes solution residual is %e\n', nlres_norm)
flowsol = xst;
%
pde=4;
it_p = 0;
%
%%% spparms('spumoni',2)
%------------------------------ nonlinear iteration
% Picard startup step
while nlres_norm>nlres0_norm*tol_nl && it_p<maxit_p,
   it_p = it_p+1;
   fprintf('\nPicard iteration number %g \n',it_p),
% compute Picard correction and update solution
%-------------------------------------------------------
%  dxns = -[Anst,Bst';Bst,sparse(np,np)]\nlres;
        % stable solve
               tt1=[ones(np1,1);zeros(np0,1)]/np1;
               tt0=[zeros(np1,1);ones(np0,1)]/np0;
        xx = -[Anst,Bst',sparse(2*nu,2); ...
              Bst,sparse(np,np),tt1,tt0; ...
              sparse(1,2*nu),tt1',zeros(1,2); ...
              sparse(1,2*nu),tt0',zeros(1,2)]\[nlres;0;0];
              dxns=xx(1:end-2); multiplier=xx(end-1:end); etoc=toc;
   xns = flowsol + dxns;
% compute residual of new solution
    N = navier_q2(xy,mv,xns);
   Anst = viscosity*A + [N, sparse(nv,nv); sparse(nv,nv), N];
   [Anst,Bst,fst,gst] = flowbc(Anst,B,f,g,xy,bound);
   np=length(gst);
   nlres = [Anst,Bst';Bst,sparse(np,np)]*xns-[fst;gst];
   nlres_norm = norm(nlres);
   nnv=length(fst); soldiff=norm(xns(1:nnv)-flowsol(1:nnv));
   fprintf('nonlinear residual is %e',nlres_norm)
   fprintf('\n    velocity change is %e',soldiff)
        p1diff=norm(xns(nnv+1:nnv+np1)-flowsol(nnv+1:nnv+np1));
        p0diff=norm(xns(nnv+np1+1:nnv+np)-flowsol(nnv+np1+1:nnv+np));
        fprintf('\n P1 pressure change is %e',p1diff)
        fprintf('\n P0 pressure change is %e\n',p0diff)
% plot flow solution
    nuv=length(fst); np1=length(xyp(:,1));
    flowplot(qmethod,xns(1:2*nu+np1),By,Bx,A,xy,xyp,x,y,bound,1,66), drawnow
    sol0=xns(2*nu+np1+1:end);
    eplot2(sol0,mv,xy,x,y,36,'pressure correction');
    pause(2)
    flowsol = xns;
% end of Picard iteration loop
    end
%
it_nl = it_p;
it_n = 0;

if nlres_norm <= nlres0_norm * tol_nl, 
   fprintf('\nfinished, nonlinear convergence test satisfied\n\n');
   nlres = nlres - [zeros(2*nv,1);(sum(nlres(2*nv+1:2*nv+np))/np)*ones(np,1)];
else
   fprintf('\nfinished, stopped on iteration counts\n\n');
end
%
%----- compute divergence error
qmethod=2; navierpost; qmethod=5;

