function x_it = m_mass(x_it,qparams)
%M_MASS mass exact matrix preconditioning operator
%   x_it = m_mass(x_it,qparams);
%   input
%          x_it       operand for preconditioning operator
%     qparams.Q       mass matrix
%   output
%          x_it       result of preconditioning operation
%
%   IFISS function: DJS; 4 August 2023.
% Copyright (c) 2010 D.J. Silvester, V. Simoncini
x_it=qparams.Q\x_it;
