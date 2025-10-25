%SOLVEX_STOKES solve two-field Stokes problem in square domain
%   IFISS scriptfile: DJS; 4 August 2023
% Copyright (c)  2023 D.J. Silvester, J. Pestana
clear variables
gohome
cd datafiles
load squarex_stokes_nobc.mat
%
%%% setup system
B=[B1;B0]; g=[g1;g0];
fprintf('imposing (enclosed flow) boundary conditions ...\n')
[Ast,Bst,fst,gst] = flowbc(A,B,f,g,xy,bound);
%
np1=length(g1); np0=length(g0); nu=length(f)/2; np=np1+np0;
%%% compute solution
tic
beta=0;
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
%-------------------------------------------------------
etoc=toc; fprintf('Stokes system solved in %8.3e seconds\n',etoc)
      fprintf('\nP1 multiplier is %8.3e',multiplier(1))
      fprintf('\nP0 multiplier is %8.3e\n',multiplier(2))
%
%%% plot solution
spc=default('uniform/nonuniform streamlines 1/2 (default uniform)',1);
flowplot(qmethod,xst(1:2*nu+np1),By,Bx,A,xy,xyp,x,y,bound,spc,35);
sol0=xst(2*nu+np1+1:end);
eplot2(sol0,mv,xy,x,y,36,'pressure correction');
fprintf('\n') 
%
%%% estimate errors
qmethod=2; stokespost; qmethod=5;

