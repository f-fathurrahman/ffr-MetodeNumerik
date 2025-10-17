function [conn,elm2edge_conn,edge2node_conn,x,y,z,exmid,eymid,ezmid, ...
    scatb_edge,scatin_edge,scatin_elm,scatb_elm,...
    huygb_elm,pmlbin_node,pmlbout_node,pml_node,pmlbout_edge] ...
    = meshgen_rectprism(delh,scat_type,lenx,leny,lenz,fsdim,pmldim,huygdim)
% Mesh generation function for "Scattering from a rectangular prism object"
% with tetrahedral elements  

% INPUT:
% delh : Element size (meter)
% scat_type : Type of scatterer, two options: 'pec' or 'diel'
% lenx : Edge length of rectangular prism along x (meter)
% leny : Edge length of rectangular prism along y (meter)
% lenz : Edge length of rectangular prism along z (meter)
% fsdim : Distance btw scatterer and inner PML boundary (meter)
% pmldim : Distance btw inner and outer PML boundaries (meter)
% huygdim : Distance btw scatterer and Huygens boundary (used if
%           scat_type='diel') (meter)
% OUTPUT:
% conn : Element-to-node connectivity matrix (size: Mx3)
% elm2edge_conn : Element-to-edge connectivity matrix (size: Mx6) 
% edge2node_conn : Edge-to-node connectivity matrix (size: Sx2)
% x,y,z  : x,y and z coordinates of all nodes in the mesh (each size: Nx1)
% exmid,eymid,ezmid : x,y and z midpoints of all edges (each size: Sx1)
% scatb_edge : Array containing the edges on the scatterer boundary
% scatin_edge : Array containing the edges inside the scatterer
% scatin_elm : Array containing the elements inside the scatterer
% scatb_elm : Array containing the elements located on the scatterer boundary
% huygb_elm : Array containing the elements located on the Huygens'
%             boundary (used if scat_type='diel')
% pmlbin_node : Array containing the nodes on the inner PML boundary
% pmlbout_node : Array containing the nodes on the outer PML boundary
% pml_node : Array containing the nodes inside the PML region
% pmlbout_edge : Array containing the edges on the outer PML boundary
%
% Note: M is the number of elements, N is the number of nodes, S is the
% number of edges
%
%
% Copyright (C) 2018, Ozlem Ozgun, Mustafa Kuzuoglu
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%    
% Contact:
% Dr. Ozlem Ozgun
% email: ozgunozlem@gmail.com
% web: http://www.ee.hacettepe.edu.tr/~ozlem/

% dimensions of cube or rectangular prism
scatmax_x = lenx/2; 
scatmax_y = leny/2;
scatmax_z = lenz/2;
scatmin_x = -lenx/2;
scatmin_y = -leny/2;
scatmin_z = -lenz/2;

huygmax_x = scatmax_x + huygdim;
huygmax_y = scatmax_y + huygdim;
huygmax_z = scatmax_z + huygdim;
huygmin_x = scatmin_x - huygdim;
huygmin_y = scatmin_y - huygdim;
huygmin_z = scatmin_z - huygdim;

pmlinmax_x = scatmax_x + fsdim;
pmlinmax_y = scatmax_y + fsdim;
pmlinmax_z = scatmax_z + fsdim;
pmlinmin_x = scatmin_x - fsdim;
pmlinmin_y = scatmin_y - fsdim;
pmlinmin_z = scatmin_z - fsdim;

pmloutmax_x = pmlinmax_x + pmldim;
pmloutmax_y = pmlinmax_y + pmldim;
pmloutmax_z = pmlinmax_z + pmldim;
pmloutmin_x = pmlinmin_x - pmldim;
pmloutmin_y = pmlinmin_y - pmldim;
pmloutmin_z = pmlinmin_z - pmldim;

% *************************************************************************
% Create coordinates
% *************************************************************************
nsx = round(abs(scatmax_x-scatmin_x)/delh)+1; %number of nodes along scat x
nsy = round(abs(scatmax_y-scatmin_y)/delh)+1; %number of nodes along scat y
nsz = round(abs(scatmax_z-scatmin_z)/delh)+1; %number of nodes along scat z

npx = round((fsdim-huygdim)/delh)+1;  % number of nodes along pml x
npy = round((fsdim-huygdim)/delh)+1;  % number of nodes along pml y
npz = round((fsdim-huygdim)/delh)+1;  % number of nodes along pml z

nhx = round(huygdim/delh)+1;  % number of nodes along huyg x
nhy = round(huygdim/delh)+1;  % number of nodes along huyg y
nhz = round(huygdim/delh)+1;  % number of nodes along huyg z

nvx = round(pmldim/delh)+1;  % number of nodes along voi x
nvy = round(pmldim/delh)+1;  % number of nodes along voi y
nvz = round(pmldim/delh)+1;  % number of nodes along voi z

tempx1 = linspace(pmloutmin_x, pmlinmin_x, nvx);
tempy1 = linspace(pmloutmin_y, pmlinmin_y, nvy);   
tempz1 = linspace(pmloutmin_z, pmlinmin_z, nvz);  

tempx2 = linspace(pmlinmin_x, huygmin_x, npx);
tempy2 = linspace(pmlinmin_y, huygmin_y, npy);   
tempz2 = linspace(pmlinmin_z, huygmin_z, npz);  

tempx3 = linspace(huygmin_x, scatmin_x, nhx);
tempy3 = linspace(huygmin_y, scatmin_y, nhy);   
tempz3 = linspace(huygmin_z, scatmin_z, nhz);  

tempx4 = linspace(scatmin_x, scatmax_x, nsx);
tempy4 = linspace(scatmin_y, scatmax_y, nsy); 
tempz4 = linspace(scatmin_z, scatmax_z, nsz);  

tempx5 = linspace(scatmax_x, huygmax_x, nhx);
tempy5 = linspace(scatmax_y, huygmax_y, nhy); 
tempz5 = linspace(scatmax_z, huygmax_z, nhz);  

tempx6 = linspace(huygmax_x, pmlinmax_x, npx);
tempy6 = linspace(huygmax_y, pmlinmax_y, npy);   
tempz6 = linspace(huygmax_z, pmlinmax_z, npz);  

tempx7 = linspace(pmlinmax_x, pmloutmax_x, nvx);
tempy7 = linspace(pmlinmax_y, pmloutmax_y, nvy);
tempz7 = linspace(pmlinmax_z, pmloutmax_z, nvz);

tempx = cat(2, tempx1, tempx2(2:end-1), tempx3, tempx4(2:end-1), ...
               tempx5, tempx6(2:end-1), tempx7);
tempy = cat(2, tempy1, tempy2(2:end-1), tempy3, tempy4(2:end-1), ...
               tempy5, tempy6(2:end-1), tempy7);
tempz = cat(2, tempz1, tempz2(2:end-1), tempz3, tempz4(2:end-1), ...
               tempz5, tempz6(2:end-1), tempz7);

[xa, ya, za] = meshgrid(tempx, tempy, tempz);

x = xa(1,:); y = ya(1,:); z = za(1,:);
 
for i = 2:size(xa, 1)
    x = cat(2, x, xa(i,:));
    y = cat(2, y, ya(i,:));
    z = cat(2, z, za(i,:));    
end 


% *************************************************************************
% Create mesh and connectivity matrix
% *************************************************************************
conn = delaunayn([x' y' z']);

N = length(x);
M = size(conn,1);

% *************************************************************************
% Find discarded elements and modify mesh if scattterer is PEC
% (make inside of scatterer empty)
% *************************************************************************
if strcmpi(scat_type, 'pec')    % PEC
    xmid = mean(x(conn),2);
    ymid = mean(y(conn),2);  
    zmid = mean(z(conn),2);    
    discard_elm = find(((xmid < (scatmax_x-1e-8)) & (xmid > (scatmin_x+1e-8))) & ...
                       ((ymid < (scatmax_y-1e-8)) & (ymid > (scatmin_y+1e-8))) & ...
                       ((zmid < (scatmax_z-1e-8)) & (zmid > (scatmin_z+1e-8))));
    
    elms0 = ones(1,M);
    elms0(discard_elm) = 0;
    nondiscard_elm = find(elms0 ~= 0);
    
    conn = conn(nondiscard_elm,:);
    
    nodes = ones(1,N);
    nodes(conn(:,1)) = 0;
    nodes(conn(:,2)) = 0;
    nodes(conn(:,3)) = 0;
    nodes(conn(:,4)) = 0;    
    nondiscard_node = find(nodes == 0);
    discard_node = find(nodes ~= 0);
    
    x = x(nondiscard_node);
    y = y(nondiscard_node);
    z = z(nondiscard_node);
    
    % renumbering
    if ~isempty(discard_node)
        connd = conn(:);
        [connd, indx] = sort(connd);
        szc = size(connd);
        connd = diff(connd);
        connd(connd>0) = 1;
        connd = cumsum(connd)+1;
        connd = [1; connd];
        connd2 = zeros(szc);
        connd2(indx) = connd;
        clear connd;
        conn = reshape(connd2, size(conn));
        clear connd2;
    end    
end

N = length(x);
M = size(conn,1);
disp(['Nelement = ' sprintf('%d', M)]);
disp(['Nnode = ' sprintf('%d', N)]);

xmid = mean(x(conn),2);
ymid = mean(y(conn),2);
zmid = mean(z(conn),2);


% *************************************************************************
% Modify nodal connectivity matrix
% *************************************************************************
[conn, ~] = modify_conn(conn, [x' y' z']);

% *************************************************************************
% Create edge connectivity matrix
% *************************************************************************
[elm2edge_conn, edge2node_conn, ~] = create_edge_conn(conn, M, N);

S = size(edge2node_conn,1);

disp(['Nedge = ' sprintf('%d', S)]);
disp(' ');


% Find the mid points of edges
exmid = (x(edge2node_conn(:,1))+x(edge2node_conn(:,2)))/2;
eymid = (y(edge2node_conn(:,1))+y(edge2node_conn(:,2)))/2;
ezmid = (z(edge2node_conn(:,1))+z(edge2node_conn(:,2)))/2;

% *************************************************************************
% Find nodes, edges and elements inside the scatterer if dielectric
% *************************************************************************
scatin_elm = [];  scatin_node = []; scatin_edge = [];
if strcmpi(scat_type, 'diel') % diel   
    scatin_elm = find(((xmid < (scatmax_x-1e-8)) & (xmid > (scatmin_x+1e-8))) & ...
                      ((ymid < (scatmax_y-1e-8)) & (ymid > (scatmin_y+1e-8))) & ...
                      ((zmid < (scatmax_z-1e-8)) & (zmid > (scatmin_z+1e-8))));
    scatin_node = find(((x < (scatmax_x+1e-8)) & (x > (scatmin_x-1e-8))) & ...
                       ((y < (scatmax_y+1e-8)) & (y > (scatmin_y-1e-8))) & ...
                       ((z < (scatmax_z+1e-8)) & (z > (scatmin_z-1e-8))));
    scatin_edge = find(((exmid < (scatmax_x+1e-8)) & (exmid > (scatmin_x-1e-8))) & ...
                       ((eymid < (scatmax_y+1e-8)) & (eymid > (scatmin_y-1e-8))) & ...
                       ((ezmid < (scatmax_z+1e-8)) & (ezmid > (scatmin_z-1e-8))));
end


% *************************************************************************
% Find scatterer boundary nodes
% *************************************************************************
scatb_node = find(((x < scatmax_x+1e-8) & (x > scatmax_x-1e-8) & ...
                  (y < scatmax_y+1e-8) & (y > scatmin_y-1e-8) & ...
                  (z < scatmax_z+1e-8) & (z > scatmin_z-1e-8)) | ...
                 ((x < scatmin_x+1e-8) & (x > scatmin_x-1e-8) & ...
                  (y < scatmax_y+1e-8) & (y > scatmin_y-1e-8) & ...
                  (z < scatmax_z+1e-8) & (z > scatmin_z-1e-8)) | ...
                 ((y < scatmax_y+1e-8) & (y > scatmax_y-1e-8) & ...
                  (x < scatmax_x+1e-8) & (x > scatmin_x-1e-8) & ...
                  (z < scatmax_z+1e-8) & (z > scatmin_z-1e-8)) | ...
                 ((y < scatmin_y+1e-8) & (y > scatmin_y-1e-8) & ...
                  (x < scatmax_x+1e-8) & (x > scatmin_x-1e-8) & ...
                  (z < scatmax_z+1e-8) & (z > scatmin_z-1e-8)) | ...
                 ((z < scatmax_z+1e-8) & (z > scatmax_z-1e-8) & ...
                  (y < scatmax_y+1e-8) & (y > scatmin_y-1e-8) & ...
                  (x < scatmax_x+1e-8) & (x > scatmin_x-1e-8)) | ...
                 ((z < scatmin_z+1e-8) & (z > scatmin_z-1e-8) & ...
                  (y < scatmax_y+1e-8) & (y > scatmin_y-1e-8) & ...
                  (x < scatmax_x+1e-8) & (x > scatmin_x-1e-8)));


% *************************************************************************
% Find Huygens boundary nodes
% *************************************************************************
huygb_node = find(((x < huygmax_x+1e-8) & (x > huygmax_x-1e-8) & ...
                   (y < huygmax_y+1e-8) & (y > huygmin_y-1e-8) & ...
                   (z < huygmax_z+1e-8) & (z > huygmin_z-1e-8)) | ...
                  ((x < huygmin_x+1e-8) & (x > huygmin_x-1e-8) & ...
                   (y < huygmax_y+1e-8) & (y > huygmin_y-1e-8) & ...
                   (z < huygmax_z+1e-8) & (z > huygmin_z-1e-8)) | ...
                  ((y < huygmax_y+1e-8) & (y > huygmax_y-1e-8) & ...
                   (x < huygmax_x+1e-8) & (x > huygmin_x-1e-8) & ...
                   (z < huygmax_z+1e-8) & (z > huygmin_z-1e-8)) | ...
                  ((y < huygmin_y+1e-8) & (y > huygmin_y-1e-8) & ...
                   (x < huygmax_x+1e-8) & (x > huygmin_x-1e-8) & ...
                   (z < huygmax_z+1e-8) & (z > huygmin_z-1e-8)) | ...
                  ((z < huygmax_z+1e-8) & (z > huygmax_z-1e-8) & ...
                   (y < huygmax_y+1e-8) & (y > huygmin_y-1e-8) & ...
                   (x < huygmax_x+1e-8) & (x > huygmin_x-1e-8)) | ...
                  ((z < huygmin_z+1e-8) & (z > huygmin_z-1e-8) & ...
                   (y < huygmax_y+1e-8) & (y > huygmin_y-1e-8) & ...
                   (x < huygmax_x+1e-8) & (x > huygmin_x-1e-8)));

% *************************************************************************
% Find PML boundary nodes
% *************************************************************************
pmlbin_node = find(((x < pmlinmax_x+1e-8) & (x > pmlinmax_x-1e-8) & ...
                   (y < pmlinmax_y+1e-8) & (y > pmlinmin_y-1e-8) & ...
                   (z < pmlinmax_z+1e-8) & (z > pmlinmin_z-1e-8)) | ...
                  ((x < pmlinmin_x+1e-8) & (x > pmlinmin_x-1e-8) & ...
                   (y < pmlinmax_y+1e-8) & (y > pmlinmin_y-1e-8) & ...
                   (z < pmlinmax_z+1e-8) & (z > pmlinmin_z-1e-8)) | ...
                  ((y < pmlinmax_y+1e-8) & (y > pmlinmax_y-1e-8) & ...
                   (x < pmlinmax_x+1e-8) & (x > pmlinmin_x-1e-8) & ...
                   (z < pmlinmax_z+1e-8) & (z > pmlinmin_z-1e-8)) | ...
                  ((y < pmlinmin_y+1e-8) & (y > pmlinmin_y-1e-8) & ...
                   (x < pmlinmax_x+1e-8) & (x > pmlinmin_x-1e-8) & ...
                   (z < pmlinmax_z+1e-8) & (z > pmlinmin_z-1e-8)) | ...
                  ((z < pmlinmax_z+1e-8) & (z > pmlinmax_z-1e-8) & ...
                   (y < pmlinmax_y+1e-8) & (y > pmlinmin_y-1e-8) & ...
                   (x < pmlinmax_x+1e-8) & (x > pmlinmin_x-1e-8)) | ...
                  ((z < pmlinmin_z+1e-8) & (z > pmlinmin_z-1e-8) & ...
                   (y < pmlinmax_y+1e-8) & (y > pmlinmin_y-1e-8) & ...
                   (x < pmlinmax_x+1e-8) & (x > pmlinmin_x-1e-8)));
               
pmlbout_node = find(((x < pmloutmax_x+1e-8) & (x > pmloutmax_x-1e-8)) | ...
                    ((x < pmloutmin_x+1e-8) & (x > pmloutmin_x-1e-8)) | ...
                    ((y < pmloutmax_y+1e-8) & (y > pmloutmax_y-1e-8)) | ...
                    ((y < pmloutmin_y+1e-8) & (y > pmloutmin_y-1e-8)) | ...
                    ((z < pmloutmax_z+1e-8) & (z > pmloutmax_z-1e-8)) | ...
                    ((z < pmloutmin_z+1e-8) & (z > pmloutmin_z-1e-8)));
               
% *************************************************************************
% Find PML inner nodes
% *************************************************************************
temp_node = find(((x < (pmlinmax_x+1e-8)) & (x > (pmlinmin_x-1e-8))) & ...
                 ((y < (pmlinmax_y+1e-8)) & (y > (pmlinmin_y-1e-8))) & ...
                 ((z < (pmlinmax_z+1e-8)) & (z > (pmlinmin_z-1e-8))));
nodes0 = ones(1,N); 
nodes0(temp_node) = 0;
pml_node = find(nodes0 == 1);


% *************************************************************************
% Find VOI and scatterer boundary edges
% *************************************************************************
nodes1 = zeros(1, N);    
nodes1(pmlbout_node) = 1;
nodes2 = zeros(1, N);    
nodes2(scatb_node) = 1;    
nodes3 = zeros(1, N);    
nodes3(huygb_node) = 1;    
             
ind1 = 1; ind2 = 1; ind3 = 1;
for i = 1:S
    n1 = edge2node_conn(i,1);
    n2 = edge2node_conn(i,2);
        
    if (nodes1(n1) == 1) && (nodes1(n2) == 1)            
        pmlbout_edge(ind1) = i;
        ind1 = ind1+1;
    end
               
    if (nodes2(n1) == 1) && (nodes2(n2) == 1)              
        scatb_edge(ind2) = i;
        ind2 = ind2+1;
    end
    
    if (nodes3(n1) == 1) && (nodes3(n2) == 1)              
        huygb_edge(ind3) = i;
        ind3 = ind3+1;
    end    
end
clear nodes1;  clear nodes2;  


% *************************************************************************
% Find scatterer and Huygens boundary elements
% *************************************************************************
indh = find(((xmid < (huygmax_x+2*delh)) & (xmid > (huygmin_x-2*delh))) & ...
            ((ymid < (huygmax_y+2*delh)) & (ymid > (huygmin_y-2*delh))) & ...
            ((zmid < (huygmax_z+2*delh)) & (zmid > (huygmin_z-2*delh))));

inds = find(((xmid < (scatmax_x)) & (xmid > (scatmin_x))) & ...
            ((ymid < (scatmax_y)) & (ymid > (scatmin_y))) & ...
            ((zmid < (scatmax_z)) & (zmid > (scatmin_z))));

elms = zeros(1,M);
elms(indh) = 1;
elms(inds) = 0;
search_elm = find(elms);

scatb_elm = []; huygb_elm = [];

ind1 = 1; ind2 = 1;
for j = 1:length(search_elm)
    elm = search_elm(j);
    
    if strcmpi(scat_type, 'pec') % PEC
        edges = zeros(1, S);    
        edges(scatb_edge) = 1;
        nodes = zeros(1, N);    
        nodes(scatb_node) = 1;

       count = 0;        
       for ie = 1:6
           if (edges(elm2edge_conn(elm,ie)) == 1)
              count = count+1;
           end
       end
        
       if (count >= 3) & ...
          ~((xmid(elm) < scatmax_x) & (xmid(elm) > scatmin_x) & ...
            (ymid(elm) < scatmax_y) & (ymid(elm) > scatmin_y) & ... 
            (zmid(elm) < scatmax_z) & (zmid(elm) > scatmin_z)),
                        
          scatb_elm(ind1,1) = elm;          
          
          nd = conn(elm,:);
          in=3;
          for i=1:4
              if (nodes(nd(i)) == 0)
                  scatb_elm(ind1,2) = nd(i);
              else 
                  scatb_elm(ind1,in) = nd(i);
                  in = in+1;
              end
          end
          if (size(scatb_elm,2) > 5) || (scatb_elm(ind1,2) == 0)
              disp('ERROR: boundary elements ...');
          end
          ind1 = ind1+1;
       end       
       
    else  % dielectric    
        edges = zeros(1, S);    
        edges(huygb_edge) = 1;
        nodes = zeros(1, N);    
        nodes(huygb_node) = 1;

       count = 0;        
       for ie = 1:6
           if (edges(elm2edge_conn(elm,ie)) == 1)
              count = count+1;
           end
       end
        
       if (count >= 3) & ...
           ~((xmid(elm) < huygmax_x+1e-8) & (xmid(elm) > huygmin_x-1e-8) & ...
             (ymid(elm) < huygmax_y+1e-8) & (ymid(elm) > huygmin_y-1e-8) & ... 
             (zmid(elm) < huygmax_z+1e-8) & (zmid(elm) > huygmin_z-1e-8))
          huygb_elm(ind2,1) = elm;          
          
          nd = conn(elm,:);
          in=3;
          for i=1:4
              if (nodes(nd(i)) == 0)
                  huygb_elm(ind2,2) = nd(i);
              else 
                  huygb_elm(ind2,in) = nd(i);
                  in = in+1;
              end
          end
          if (size(huygb_elm,2) > 5) || (huygb_elm(ind2,2) == 0)
              disp('ERROR: boundary elements ...');
          end
          ind2 = ind2+1;
       end 
              
    end  
end


% *************************************************************************
% Plot results
% *************************************************************************
plot_flag = 1; % flag showing whether the mesh will be plotted or not
if plot_flag
%     figure; hold on
%     tetramesh(conn,[x' y' z']); camorbit(20,0);
%     axis equal tight;
%     set(gcf,'Color',[1 1 1]);
%     xlabel('x (m)');  ylabel('y (m)'); zlabel('y (m)');

    figure; hold on
    plot3(exmid(scatb_edge),eymid(scatb_edge),ezmid(scatb_edge),'b.')
    plot3(x(pmlbin_node),y(pmlbin_node),z(pmlbin_node),'m.')
    plot3(x(pmlbout_node),y(pmlbout_node),z(pmlbout_node),'k.') 
    plot3(x(pml_node),y(pml_node),z(pml_node),'y.')
    plot3(xmid(scatin_elm),ymid(scatin_elm),zmid(scatin_elm),'c.')
    if ~isempty(scatb_elm), plot3(xmid(scatb_elm(:,1)),ymid(scatb_elm(:,1)),zmid(scatb_elm(:,1)),'g.'); end
    if ~isempty(huygb_elm), plot3(xmid(huygb_elm(:,1)),ymid(huygb_elm(:,1)),zmid(huygb_elm(:,1)),'g.'); end   
    axis equal tight;
    set(gcf,'Color',[1 1 1]);
    xlabel('x (m)');  ylabel('y (m)'); zlabel('y (m)');
    view(-30,30)
end
  
% keyboard

%%
% *************************************************************************
% *************************************************************************
% *************************************************************************

function [elm2edge_conn, edge2node_conn, node2element] = create_edge_conn(conn, M, N)
% It creates edge connectivity arrays from element-to-node conectivity array,
% for tetrahedral elements. 
% INPUT: 
% conn: node IDs for each element (nodal connectivity matrix) (Size: M x 4)
% OUTPUT: 
% elm2edge_conn: edge IDs for each element (Size: M x 6)
% edge2node_conn: node IDs for each edge (Size: Nedge x 2)
% node2element: elements connected to each node (Size: N x something)

TETRA_NUM_OF_NODES = 4;
TETRA_NUM_OF_EDGES = 6;

local_edge_def = [1 2; ... % edgeID=1
                  1 3; ... % edgeID=2
                  1 4; ... % edgeID=3
                  2 3; ... % edgeID=4
                  2 4; ... % edgeID=5
                  3 4];    % edgeID=6
% temporary, order of each row may be reversed

node2element = struct('elm', []);
elm2edge_conn = zeros(M, TETRA_NUM_OF_EDGES);

% Form node2element array (elements connected to each node)
for i = 1:N
    [rowind, ~] = find(conn == i);
    node2element(i).elm = rowind;
    node2element(i).elm = sort(node2element(i).elm);
end

% Form elm2edge_conn & edge2node_conn arrays
ecounter = 1; % edge counter
for i = 1:M
    for j = 1:TETRA_NUM_OF_EDGES
        if (elm2edge_conn(i,j) == 0)
            elm2edge_conn(i,j) = ecounter;
            
            % Find the two local nodes from local edge ID
            localn1 = local_edge_def(j, 1);
            localn2 = local_edge_def(j, 2);
                        
            globn1 = conn(i, localn1);
            globn2 = conn(i, localn2);
            
            % Form edge2node_conn array
            edge2node_conn(ecounter, 1) = min(globn1, globn2);
            edge2node_conn(ecounter, 2) = max(globn1, globn2);
            
            for k = 1:length(node2element(globn1).elm)
                if ~isempty(find(node2element(globn2).elm == node2element(globn1).elm(k))) & (node2element(globn1).elm(k) ~= i)
                   newi = node2element(globn1).elm(k);

                   % Find the local nodes from global nodes for element elm=newi
                   localn1 = find(conn(newi,:) == globn1);
                   localn2 = find(conn(newi,:) == globn2);
	
                   newj = LocalEdge_From_Local_Nodes(localn1, localn2);
                   elm2edge_conn(newi,newj) = ecounter;                   
                end               
            end
            ecounter = ecounter+1;
        end
    end
end


%%
function [ledgeID] = LocalEdge_From_Local_Nodes(node1, node2)
% It returns the local edge ID from two local nodes for tetrahedra
	
if (((node1==1) & (node2==2)) | ((node1==2) & (node2==1)))
	ledgeID = 1;
elseif (((node1==1) & (node2==3)) | ((node1==3) & (node2==1)))
	ledgeID = 2;
elseif (((node1==1) & (node2==4)) | ((node1==4) & (node2==1)))
	ledgeID = 3;
elseif (((node1==2) & (node2==3)) | ((node1==3) & (node2==2)))
	ledgeID = 4;
elseif (((node1==4) & (node2==2)) | ((node1==2) & (node2==4)))
	ledgeID = 5;
elseif (((node1==3) & (node2==4)) | ((node1==4) & (node2==3)))
	ledgeID = 6;
end

%%
function [conn, volume] = modify_conn(conn, coords)
% This function modifies connectivity matrix such that the volume of 
% tetrahedron is nonnegative, for each element

M = size(conn,1);
P = perms([4 3 2 1]);  % permutation array
volume = zeros(1,M);

for elm = 1:M
    for ip = 1:size(P,1)
        
        n1 = conn(elm,P(ip,1));
        n2 = conn(elm,P(ip,2));
        n3 = conn(elm,P(ip,3));
        n4 = conn(elm,P(ip,4));
        
        volume(elm) = tetrahedron_volume(coords(n1,:), coords(n2,:), coords(n3,:), coords(n4,:));

        if volume(elm) > 0 
           conn(elm,:) = [n1 n2 n3 n4];
           break;
        end
    end         
end

%%
function [vol] = tetrahedron_volume(a, b, c, d)
% This function calculates the signed volume of a tetrahedron
% whose vertices are a,b,c,d

ax = a(1); ay = a(2); az = a(3);
bx = b(1); by = b(2); bz = b(3);
cx = c(1); cy = c(2); cz = c(3);
dx = d(1); dy = d(2); dz = d(3);

vol = det([1 1 1 1; ax bx cx dx; ay by cy dy; az bz cz dz])/6;
