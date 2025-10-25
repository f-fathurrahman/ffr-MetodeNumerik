%SOLVEX_STEP_STOKES solve two-field Stokes problem in step domain
%   IFISS scriptfile: DJS;  4 August 2023
% Copyright (c) 2023 D.J. Silvester, J. Pestana
clear variables
gohome
cd datafiles
load stepx_stokes_nobc.mat
%
%%% setup system
B=[B1;B0]; g=[g1;g0];
fprintf('imposing essential boundary conditions ...\n')
[Ast,Bst,fst,gst] = flowbc(A,B,f,g,xy,bound);
%
np1=length(g1); np0=length(g0); nu=length(f)/2; np=np1+np0;
%%% compute solution
tic
beta=0;
%-------------------------------------------------------
% xref=[Ast,Bst';Bst,sparse(np,np)]\[fst;gst];
% stable solve
      tt=[ones(np1,1);-ones(np0,1)]/np;
xx = [Ast,Bst',sparse(2*nu,1); ...
      Bst,sparse(np,np),tt; ...
      sparse(1,2*nu),tt',zeros(1,1)]\[fst;gst;0];
      xst=xx(1:end-1); multiplier=xx(end); etoc=toc;
%-------------------------------------------------------
etoc=toc; fprintf('Stokes system solved in %8.3e seconds\n',etoc)
      fprintf('P1P0 multiplier is %8.3e\n',multiplier)

%
%%% plot solution
flowplotl(qmethod,xst(1:2*nu+np1),By,Bx,A,xy,xyp,x,y,bound,35);
sol0=xst(2*nu+np1+1:end);
eplotl(sol0,mv,xy,x,y,36,'pressure correction');
fprintf('\n') 
%
%%% estimate errors
%qmethod=2; stokespost; qmethod=5;

