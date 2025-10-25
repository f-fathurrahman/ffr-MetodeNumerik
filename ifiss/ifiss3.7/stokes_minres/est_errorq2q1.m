function errorest = est_errorq2q1(xmr,eparams)
%EST_ERRORQ2Q1 computes energy error estimate for Q2-Q1 solution 
%   errorest = est_errorq2q1(xmr,eparams);
%   input
%          xmr        Q2-Q1 solution iterate
%      eparams        structure for error estimator:
%          .ae        elementwise Poisson problem matrices
%          .xy        vertex coordinate vector  
%          .mv        element mapping matrix
%      .mbound        element edge boundary matrix 
%   output
%     errorest        energy error estimate
%
%   IFISS function: DJS; 21 May 2010, 8 August 2023
% Copyright (c) 2010 D.J. Silvester, Qifeng Liao

%% unpack eparams structure
ae=eparams.ae; xy=eparams.xy; mv=eparams.mv; mbound=eparams.mbound;
neumannb=eparams.neumannb; eex=eparams.eex; hx=eparams.hx; hy=eparams.hy;
mp=eparams.mp;
[error_x,error_y] = algpost_q2q1(xmr,ae,xy,mv,mp,mbound,neumannb,eex,hx,hy);
error_div = q2div(xy,mv,xmr,0);
errorest=norm(sqrt(error_x+error_y+error_div.^2),2);
