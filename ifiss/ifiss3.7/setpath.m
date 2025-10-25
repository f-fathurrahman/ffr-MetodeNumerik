%SETPATH sets IFISS search path 
%   IFISS scriptfile: DJS; 29 December 2024.
% Copyright (c) 2009 D.J. Silvester, H.C. Elman, A. Ramage 

% initilise matlab path for ./ifiss/ directories
gohome, 
warning('off', 'all')
addpath(genpath(pwd),'-end')
%fprintf('Resolving these path conflicts ... done')
warning('on','all')
%

% fix iterapp bug (affects pcg & minres)
if strncmp(version,'7.0.4',5), 
error('Oops ... unsupported version of MATLAB')
end

% include old version of ilu 
if strncmp(version,'6.5',3) | strncmp(version,'7.0',3) | strncmp(version,'7.1.',4) | ...
    strncmp(version,'7.2',3) | strncmp(version,'7.3',3), 
error('Oops ... unsupported version of MATLAB')
end

% fix MATLAB UMFPACK bug
if strncmp(version,'6.5',3) | strncmp(version,'7.0',3), 
spparms('default'),spparms('piv_tol',1), end

