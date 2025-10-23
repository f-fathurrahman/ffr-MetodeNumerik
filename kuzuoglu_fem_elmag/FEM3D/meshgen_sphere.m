function [conn,elm2edge_conn,edge2node_conn,x,y,z,exmid,eymid,ezmid,...
    scatb_edge,scatin_edge,scatin_elm,scatb_elm,...
    huygb_elm,pmlbin_node,pmlbout_node,pml_node,pmlbout_edge] ...
    = meshgen_sphere(delh,scat_type,rs,fsdim,pmldim,huygdim)
% Mesh generation function for "Scattering from a spherical object"
% with tetrahedral elements  

% INPUT:
% delh : Element size (meter)
% scat_type : Type of scatterer, two options: 'pec' or 'diel'
% rs : Radius of sphere (meter)
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

rpmlin = rs+fsdim;       % radius of inner PML boundary
rpmlout = rpmlin+pmldim; % radius of outer PML boundary
rhuyg = rs+huygdim;      % radius of Huygens boundary

% fit to delh incrementation
rs = round(rs/delh)*delh; 
rpmlout = round(rpmlout/delh)*delh; 
rpmlin = round(rpmlin/delh)*delh; 
rhuyg = round(rhuyg/delh)*delh; 
             
% *************************************************************************
% Create coordinates layer by layer
% *************************************************************************
x = []; y = []; z = [];

if strcmpi(scat_type, 'diel') % dielectric
    x = 0; y = 0; z= 0; % add center to mesh
end

for r = delh:delh:rpmlout
    if (strcmpi(scat_type, 'diel')) || ((strcmpi(scat_type, 'pec')) && (r >= rs-1e-8))
        % create a single sphere with radius r
        x1 = []; y1 = []; z1 = [];
        zi = 0;
        firstflag = 1;
        zflag = 1;
        while (zi < r) & zflag
            if firstflag
                zi = 0;
            else
                ziold = zi;
                alfas = 2*asin(delh*0.5/radiusa);
                zi = zi+radiusa*tan(alfas);
            end
            
            if (zi >= r)
                zi = ziold + 2*(zi-ziold)/3;
                if (zi >= r)
                    zi = ziold + 2*(r-ziold)/3;
                end
                zflag = 0;
            end
            
            radiusa = sqrt(r^2-zi^2);
            alfac = 2*asin(delh*0.5/radiusa);
            
            N = floor(pi/alfac);
            alfa = cat(2, (pi/N)*(0:N), -fliplr((pi/N)*(1:N-1)));
            
            x0 = radiusa*cos(alfa);
            y0 = radiusa*sin(alfa);
            z0 = zi*ones(size(x0));
            
            x1 = cat(2, x1, x0);
            y1 = cat(2, y1, y0);
            z1 = cat(2, z1, z0);
            
            if ~firstflag
                x1 = cat(2, x1, x0);
                y1 = cat(2, y1, y0);
                z1 = cat(2, z1, -z0);
            end
            firstflag = 0;
        end % while
        
        x1 = cat(2, x1, 0, 0);
        y1 = cat(2, y1, 0, 0);
        z1 = cat(2, z1, r, -r);
        
        x = cat(2, x, x1);
        y = cat(2, y, y1);
        z = cat(2, z, z1);

        if (r < rs+1e-8) && (r > rs-1e-8)
            scatb_x = x1; scatb_y = y1; scatb_z = z1;
        elseif (r < rhuyg+1e-8) && (r > rhuyg-1e-8)
            huygb_x = x1; huygb_y = y1; huygb_z = z1;
        elseif (r < rpmlin+1e-8) && (r > rpmlin-1e-8)
            pmlbin_x = x1; pmlbin_y = y1; pmlbin_z = z1;
        elseif (r < rpmlout+1e-8) && (r > rpmlout-1e-8)
            pmlbout_x = x1; pmlbout_y = y1; pmlbout_z = z1;
        end
    end
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
    discard_elm = find(sqrt(xmid.^2 + ymid.^2 + zmid.^2) < rs-1e-6);    
    
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
[conn, volume] = modify_conn(conn, [x' y' z']);

% *************************************************************************
% Create edge connectivity matrix
% *************************************************************************
[elm2edge_conn, edge2node_conn, ~] = create_edge_conn(conn, M, N);

Nedge = size(edge2node_conn,1);

disp(['Nedge = ' sprintf('%d', Nedge)]);
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
    scatin_elm = find(sqrt(xmid.^2 + ymid.^2 + zmid.^2) < rs);
    scatin_node = find(sqrt(x.^2 + y.^2 + z.^2) < rs+1e-6);
    scatin_edge = find(sqrt(exmid.^2 + eymid.^2 + ezmid.^2) < rs+1e-6);
end

% *************************************************************************
% Find scatterer boundary nodes
% *************************************************************************
[scatb_node,~] = knnsearch([x' y' z'],[scatb_x' scatb_y' scatb_z'],'k',1);
scatb_node = sort(unique(scatb_node));

% *************************************************************************
% Find Huygens boundary nodes
% *************************************************************************
[huygb_node,~] = knnsearch([x' y' z'],[huygb_x' huygb_y' huygb_z'],'k',1);
huygb_node = sort(unique(huygb_node));

% *************************************************************************
% Find PML boundary nodes
% *************************************************************************
[pmlbin_node,~] = knnsearch([x' y' z'],[pmlbin_x' pmlbin_y' pmlbin_z'],'k',1);
pmlbin_node = sort(unique(pmlbin_node));

[pmlbout_node,~] = knnsearch([x' y' z'],[pmlbout_x' pmlbout_y' pmlbout_z'],'k',1);
pmlbout_node = sort(unique(pmlbout_node));

% *************************************************************************
% Find PML inner nodes
% *************************************************************************
temp_node = find(sqrt(x.^2 + y.^2 + z.^2) < rpmlin+1e-6);
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

d = sqrt(exmid.^2 + eymid.^2 + ezmid.^2);
search_edge = find(((d < rpmlout+1*delh) & (d > rpmlout-1*delh)) | ...
                   ((d < rs+1*delh) & (d > rs-1*delh)) | ...
                   ((d < rhuyg+1*delh) & (d > rhuyg-1*delh)));

ind1 = 1; ind2 = 1; ind3 = 1;
for ii = 1:length(search_edge)
    i = search_edge(ii);
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
d = sqrt(xmid.^2 + ymid.^2 + zmid.^2);
search_elm = find(((d < rhuyg+2*delh) & (d > rhuyg-2*delh)) | ...
                  ((d < rs+2*delh) & (d > rs-2*delh)));

scatb_elm = []; huygb_elm = [];

ind1 = 1; ind2 = 1;
for j = 1:length(search_elm)
    elm = search_elm(j);
    
    if strcmpi(scat_type, 'pec') % PEC
        edges = zeros(1, Nedge);    
        edges(scatb_edge) = 1;
        nodes = zeros(1, N);    
        nodes(scatb_node) = 1;

       count = 0;        
       for ie = 1:6
           if (edges(elm2edge_conn(elm,ie)) == 1)
              count = count+1;
           end
       end
        
       if (count >= 3) && (sqrt(xmid(elm)^2 + ymid(elm)^2 + zmid(elm)^2) > rs)
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
        edges = zeros(1, Nedge);    
        edges(huygb_edge) = 1;
        nodes = zeros(1, N);    
        nodes(huygb_node) = 1;

       count = 0;        
       for ie = 1:6
           if (edges(elm2edge_conn(elm,ie)) == 1)
              count = count+1;
           end
       end
        
       if (count >= 3) && (sqrt(xmid(elm)^2 + ymid(elm)^2 + zmid(elm)^2) > rhuyg)
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
