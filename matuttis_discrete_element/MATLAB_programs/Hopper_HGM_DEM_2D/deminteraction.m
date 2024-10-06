function [acx,acy,acf] =deminteraction(...
  col_list,col_list_len,dt,particleX,particleY,...
  particleSide,PX,PY,VX,VY,VF,freeParticle,m,mi)
% PUPOSE: driver to call the computation of the 
%  overlap geometry and from that the 
%  interaction forces
% CAVEAT: This is the most CPU-time-consuming routine 
%  for dense systems; speeding up other functions will
%  hardly have any impact
% LITERATURE: sec. 7.2,7.3 
% REVISION HISTORY: Shi Han Ng, HG Matuttis 15-May-2014
%
global gravity
global oldRealList
global oldRealListLen
global oldTtanfList
global oldAreaList
global oldVreltList

particleNum=length(particleSide);
acx=zeros(particleNum,1);
acy=zeros(particleNum,1);
acf=zeros(particleNum,1);

realList    = zeros(2,8*particleNum);
realListLen = 0;
ttanfList   = zeros(1,8*particleNum);
areaList    = zeros(1,8*particleNum);
vreltList   = zeros(1,8*particleNum);

% Loop over the collision list which
% contains the particle pairs, i.e. the
% respective two indices of each of the interacting
% pairs
for i = 1:col_list_len
  % Reorder the indices of the interacting pairs 
  % so that the index for the first particle
  % is smaller than for the second particle 
  if col_list(1,i) < col_list(2,i)
    polyA = col_list(1,i);
    polyB = col_list(2,i);
  else
    polyA = col_list(2,i);
    polyB = col_list(1,i);
  end
 
% If-condition serves to skip physically irrelevant 
% computation of overlaps between walls which may
% be geometrically problematic
  if (freeParticle(polyA)|freeParticle(polyB))
    [overlap_geometry] = ...
      overlap(particleX(1:particleSide(polyA),polyA),...
              particleY(1:particleSide(polyA),polyA),...
              PX(polyA),PY(polyA),...
              particleX(1:particleSide(polyB),polyB),...
              particleY(1:particleSide(polyB),polyB),...
              PX(polyB),PY(polyB)); %#ok<NASGU>

    % FORCE COMPUTATION
    if( overlap_geometry(1) > 0 )
      % Add the pair into the collision list.
      realListLen = realListLen + 1;
      realList(:,realListLen) = [polyA;polyB];

      % Check the pair if they are in the old list.
      oldTtanf = 0;
      oldArea  = 0;
      oldVrelt = 0;
      for iOldReal = 1:oldRealListLen;
        if ( (oldRealList(1,iOldReal) == polyA) && ...
             (oldRealList(2,iOldReal) == polyB) )
           oldTtanf = oldTtanfList(iOldReal);
           oldArea  = oldAreaList(iOldReal);
           oldVrelt = oldVreltList(iOldReal);
           break;
        end
      end

      [resulting_forces,ttanfList(realListLen),areaList(realListLen),...
      vreltList(realListLen)] = ...
      findForce(overlap_geometry,...
      VX(polyA),VY(polyA),VF(polyA),m(polyA),...
      PX(polyA),PY(polyA),mi(polyA),...
      VX(polyB),VY(polyB),VF(polyB),m(polyB),...
      PX(polyB),PY(polyB),mi(polyB),...
      oldTtanf,oldArea,oldVrelt,dt);

      acx(polyA)=acx(polyA)+resulting_forces(1);
      acy(polyA)=acy(polyA)+resulting_forces(2);
      acf(polyA)=acf(polyA)+resulting_forces(3);

      acx(polyB)=acx(polyB)-resulting_forces(1);
      acy(polyB)=acy(polyB)-resulting_forces(2);
      acf(polyB)=acf(polyB)+resulting_forces(5);
    end
  end
end

% compute accelerations from forces
acx=acx./m(:);
acy=acy./m(:);
acf=acf./mi(:);
% add gravitational acceleration
acy=acy-gravity; 
% the accelerations for the walls will be set to zero
% with the use of the freeParticle-array in gearCorrector
% so that the walls won't move

oldRealList    = realList;
oldRealListLen = realListLen;
oldTtanfList   = ttanfList;
oldAreaList    = areaList;
oldVreltList   = vreltList;

return
