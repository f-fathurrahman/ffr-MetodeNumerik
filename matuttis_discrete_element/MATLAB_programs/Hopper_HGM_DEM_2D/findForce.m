function[resulting_forces,old_ttanf,old_area,old_vrelt] =...
  findForce(overlap_geometry,vx0,vy0,vf0,mass0,pcx0,...
  pcy0,mi0,vx1,vy1,vf1,mass1,pcx1,pcy1,mi1,old_ttanf,old_area,old_vrelt,dt)
% PURPOSE: Computing forces between two particles.
% USAGE: Call with the input computed from overlap,
%   pass the results to the corrector
% Literature: section 7.3
%  Detailed Numbers for formulae and pages 
%  are given in the listing
% REVISION HISTORY: Shi Han Ng, HG Matuttis 15-May-2014
%

global emodul
global mu
global gamma

% Extract input-data from input-vector
area=overlap_geometry(1);
% force point
scx =overlap_geometry(2); 
scy =overlap_geometry(3);
% normal, simplified version of Nr. 8, p 250
% should not be normalized, but the length
% contains the contact length
nx  =overlap_geometry(4);
ny  =overlap_geometry(5);

contactlength=sqrt(nx*nx+ny*ny);
if (abs(contactlength-1.0)<1e-8)
  disp('Probably wrong version of the overlap computation.')
  disp('Please verify!!!')
  return
end

% nx, ny are not normalized to 1, but 
% contain the length of the contact line
if (contactlength~=0)
  lambda=1/contactlength;
  nx=nx*lambda;
  ny=ny*lambda;
end

% vectors from centers of mass (pcx0,pcy0)
% and (pcx1,pcy1) to the contact point (scx,scy)
r0x=scx-pcx0;
r0y=scy-pcy0;
r1x=scx-pcx1;
r1y=scy-pcy1;

r0abs=sqrt(r0x^2+r0y^2);
r1abs=sqrt(r1x^2+r1y^2);
% Inverse of eq. (7.6)
icl =(r0abs+r1abs)/(4*r0abs*r1abs);

% Elastic contact force (7.7)
fn_rep=emodul*area*icl;
% If a cohesive force is implemented (sec. 7.3.5), 
% it is advisable be in multiples of the Young's 
% modulus % so that one can estimate the strength
% fn_coh=-0.02 *emodul* contactlength; 
fn_coh=0.0; % currently, no cohesion implemented

% relative velocities, Fig. 7.17
vrelx=(vx0-vx1) - (-r1y*vf1 + r0y*vf0);
vrely=(vy0-vy1) - ( r1x*vf1 - r0x*vf0);
% tangential relative velocities
vrelt=-vrelx*ny + vrely*nx;
% normal relative velocities may be needed for
% alternative definition of the damping according 
% to  (7.11) 
%vreln= vrelx*nx + vrely*ny;

% effective mass/ reduced mass of the particle pair
meff    =(mass0*mass1)/(mass0+mass1);
% and tangential reduced mass mass, see p. 144 
meff_tan=1/(1/mass0+1/mass1+r0abs^2/mi0+r1abs^2/mi1);

% Normal damping: 7.3.4
damp=gamma*sqrt(meff*emodul);

% Damping in normal direction (7.13)
fni=-damp*(old_area-area)*icl/dt;
%Alternatively use  (7.11) 
% fni=-damp*vreln

% Avoid jumps of the force at approach
% due to large changes in the damping force
% as discussed on p. 227
if (fni>2*fn_rep)
  fni=abs(fn_rep)*sign(fni);
end

% Realistic damping: the cohesive part during
% separation should be cut off as discussed
% on p. 277
if (fn_rep*(fn_rep+fni)<0)
  if (abs(fni)>abs(fn_rep))
    fni=-fn_rep;
  end
end
fn=fni+fn_rep+fn_coh;

% Upper limit Coulomb friction due to 
% normal force and friction coefficient
ft=abs(mu*fn);

% Cundall-Strack model 3.4.1
% Tangential "spring constant", exact for
% two-particle contacts of spheres, larger 
% constants can also be used, but they may 
% necessitate smaller timesteps:
k_tangential=emodul*2.0/7.0;

dttanf=-k_tangential*vrelt*dt;
% If there was no normal force in the
% previous timestep, there should also
% not be a tangential force
if (old_area==0)
  if (old_ttanf~=0)
    disp('Bug in lists!!!');
    keyboard
  end
end

% Cundall-Strack friction:
% Incrementation of the tangential force
ttanf=old_ttanf+dttanf;
% Cut-off for Coulomb-friction for the hysteretic force
if (abs(ttanf)>ft)
  ttanf=abs(ft)*sign(ttanf);
end

% Numerical prefactor for the tangential Damping 
pref_y_v_tan=0.25*sqrt(abs(k_tangential*meff_tan)); 

% add tangential damping proportional to -V.
acting_ttanf=ttanf - pref_y_v_tan*vrelt;  
% Cut-off for Coulomb-friction for the sum
% of hysteretic force and tangential damping
if (abs(acting_ttanf)>ft)
  acting_ttanf=abs(ft)*sign(acting_ttanf);
end
%(Only one cut-off is mentioned in 3.4.1)

% Forces with components in x- and y-direction
%(prefactors -,+ are a result of the
% 'two-dimensional' vector product)
fxt0=fn*nx-acting_ttanf*ny;
fyt0=fn*ny+acting_ttanf*nx;

% Torques p. 240 (different for the two particles
% if the distances of the contact point to the
% centers of mass are different
fft0= (r0x*fyt0-r0y*fxt0);
fft1=-(r1x*fyt0-r1y*fxt0);

% save parameters for use in the next timestep
old_area =area;
old_vrelt=vrelt;
old_ttanf=ttanf;

% copy forces and torques into the vector
% array for use in the current timestep
resulting_forces(1)=fxt0;
resulting_forces(2)=fyt0;
resulting_forces(3)=fft0;
resulting_forces(5)=fft1;

return

