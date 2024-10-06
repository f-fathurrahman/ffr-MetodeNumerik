function [particleX,particleY]=updatevertices(particleSide,...
  PXOld,PYOld,PFOld,PX,PY,PF,particleX,particleY)
% PURPOSE: Rotation with respect to the center of mass
%  and translation according to the difference
%  between the current and the previous timestep
% USAGE: If a predictor-corrector integrator is used,
%  in timestep n it must be called after the predictor-step, 
%  as the corrected coordinates of timestep (n-1) and the
%  predicted coordinates of timestep n are used.
% ALGORITHM:
%  1. Shift the center of mass into the origin,
%  2. Rotate the particle with an angle which is
%     the difference between the old and the new orientation
%  3. Translate the corners to the new center
%     of mass (predicted coordinate)
% REVISION HISTORY: Shi Han Ng, HG Matuttis 15-May-2014

particleNum=length(PX);

particleRelX=zeros(size(particleX));
particleRelY=zeros(size(particleY));

for i=1:particleNum
% Relative Position=Center of Mass - Absolute Position
  particleRelX(1:particleSide(i),i) = PXOld(i) - particleX(1:particleSide(i),i);
  particleRelY(1:particleSide(i),i) = PYOld(i) - particleY(1:particleSide(i),i);
  ps=particleSide(i);
  % ROTATION AND PARTICLES UPDATES
  [particleX(1:ps,i),particleY(1:ps,i)]=...
  rotation(particleRelX(1:ps,i),particleRelY(1:ps,i),...
  PX(i),PY(i),PF(i)-PFOld(i));
end

end

