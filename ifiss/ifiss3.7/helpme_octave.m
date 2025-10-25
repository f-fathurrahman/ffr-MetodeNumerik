%HELPME_OCTAVE  Octave specific issues 
%   IFISS scriptfile: DJS; 29 December 2024.
% Copyright (c) 2012 D.J. Silvester, H.C. Elman, A. Ramage
fprintf(' \n');
fprintf(' IFISS 3.7 \n')
fprintf(' The Octave package is identical to the MATLAB package\n');
fprintf(' \n');
fprintf(' Here is a list of minor "issues" with testing with Octave 9.3 \n');
fprintf(' *  There are a host of warning messages associated with \n')
fprintf('    incompatibility (mostly due to sloppy syntax) ... \n');
fprintf('    typing "warning off all" effectively supresses these \n');
fprintf(' *  The generation of movies requires additional Octave toolbox\n');
fprintf('    functionality (there is no "WriteVideo" in standard Octave)\n');
fprintf(' \n');
fprintf(' Please report any other compatibility problems that you encounter\n') 
fprintf(' to the Octave Bug Report website.\n\n');
