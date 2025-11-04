%==========================================================================
% Display cell view of problem geometry
figure_number = figure;

for ind =1:size(material_types,2)
    matcol(ind,:) = material_types(ind).color;
end

m3d  = ones(nx+2, ny+2, nz+2);
tm3d = zeros(nx, ny, nz);
m3d(2:nx+1, 2:ny+1, 2:nz+1) = material_3d_space;
sz = [nx ny nz];

tm3d = m3d(2:nx+1, 2:ny+1, 2:nz+1) ...
    - m3d(1:nx, 2:ny+1, 2:nz+1);
IND = find(tm3d>0);
[i1,j1,k1] = ind2sub(sz,IND);

tm3d = m3d(2:nx+1, 2:ny+1, 2:nz+1) ...
    - m3d(3:nx+2, 2:ny+1, 2:nz+1);
IND = find(tm3d>0);
[i2,j2,k2] = ind2sub(sz,IND);

tm3d = m3d(2:nx+1, 2:ny+1, 2:nz+1) ...
    - m3d(2:nx+1, 1:ny, 2:nz+1);
IND = find(tm3d>0);
[i3,j3,k3] = ind2sub(sz,IND);

tm3d = m3d(2:nx+1, 2:ny+1, 2:nz+1) ...
    - m3d(2:nx+1, 3:ny+2, 2:nz+1);
IND = find(tm3d>0);
[i4,j4,k4] = ind2sub(sz,IND);

tm3d = m3d(2:nx+1, 2:ny+1, 2:nz+1) ...
    - m3d(2:nx+1, 2:ny+1, 1:nz);
IND = find(tm3d>0);
[i5,j5,k5] = ind2sub(sz,IND);

tm3d = m3d(2:nx+1, 2:ny+1, 2:nz+1) ...
    - m3d(2:nx+1, 2:ny+1, 3:nz+2);
IND = find(tm3d>0);
[i6,j6,k6] = ind2sub(sz,IND);
clear tm3d;

xx = node_coordinates_xe;
yy = node_coordinates_ye;
zz = node_coordinates_ze;
siz = [nx+2 ny+2 nz+2];

len = size(i1,1);
vx = [xx(i1).' xx(i1).' xx(i1).' xx(i1).'];
vy = [yy(j1).' yy(j1+1).'  yy(j1+1).' yy(j1).'];
vz = [zz(k1).' zz(k1).' zz(k1+1).' zz(k1+1).'];
vx  = reshape(vx.',4*len,1);    
vy  = reshape(vy.',4*len,1);    
vz  = reshape(vz.',4*len,1);    
v1 = [vx vy vz];
f1 = [(1:4:4*(len-1)+1).' (2:4:4*(len-1)+2).' (3:4:4*(len-1)+3).' (4:4:4*(len-1)+4).'];

IND = sub2ind(siz,i1+1, j1+1, k1+1);
RGB1 = matcol(m3d(IND),:);

len = size(i2,1);
vx = [xx(i2+1).' xx(i2+1).' xx(i2+1).' xx(i2+1).'];
vy = [yy(j2).' yy(j2+1).'  yy(j2+1).' yy(j2).'];
vz = [zz(k2).' zz(k2).' zz(k2+1).' zz(k2+1).'];
vx  = reshape(vx.',4*len,1);    
vy  = reshape(vy.',4*len,1);    
vz  = reshape(vz.',4*len,1);    
v2 = [vx vy vz];
f2 = [(1:4:4*(len-1)+1).' (2:4:4*(len-1)+2).' (3:4:4*(len-1)+3).' (4:4:4*(len-1)+4).'];

IND = sub2ind(siz,i2+1, j2+1, k2+1);
RGB2 = matcol(m3d(IND),:);

len = size(i3,1);
vx = [xx(i3).' xx(i3+1).' xx(i3+1).' xx(i3).'];
vy = [yy(j3).' yy(j3).' yy(j3).' yy(j3).'];
vz = [zz(k3).' zz(k3).' zz(k3+1).' zz(k3+1).'];
vx  = reshape(vx.',4*len,1);    
vy  = reshape(vy.',4*len,1);    
vz  = reshape(vz.',4*len,1);    
v3 = [vx vy vz];
f3 = [(1:4:4*(len-1)+1).' (2:4:4*(len-1)+2).' (3:4:4*(len-1)+3).' (4:4:4*(len-1)+4).'];

IND = sub2ind(siz,i3+1, j3+1, k3+1);
RGB3 = matcol(m3d(IND),:);

len = size(i4,1);
vx = [xx(i4).' xx(i4+1).' xx(i4+1).' xx(i4).'];
vy = [yy(j4+1).' yy(j4+1).' yy(j4+1).' yy(j4+1).'];
vz = [zz(k4).' zz(k4).' zz(k4+1).' zz(k4+1).'];
vx  = reshape(vx.',4*len,1);    
vy  = reshape(vy.',4*len,1);    
vz  = reshape(vz.',4*len,1);    
v4 = [vx vy vz];
f4 = [(1:4:4*(len-1)+1).' (2:4:4*(len-1)+2).' (3:4:4*(len-1)+3).' (4:4:4*(len-1)+4).'];

IND = sub2ind(siz,i4+1, j4+1, k4+1);
RGB4 = matcol(m3d(IND),:);

len = size(i5,1);
vx = [xx(i5).' xx(i5).' xx(i5+1).' xx(i5+1).'];
vy = [yy(j5).' yy(j5+1).' yy(j5+1).' yy(j5).'];
vz = [zz(k5).' zz(k5).' zz(k5).' zz(k5).'];
vx  = reshape(vx.',4*len,1);    
vy  = reshape(vy.',4*len,1);    
vz  = reshape(vz.',4*len,1);    
v5 = [vx vy vz];
f5 = [(1:4:4*(len-1)+1).' (2:4:4*(len-1)+2).' (3:4:4*(len-1)+3).' (4:4:4*(len-1)+4).'];

IND = sub2ind(siz,i5+1, j5+1, k5+1);
RGB5 = matcol(m3d(IND),:);

len = size(i6,1);
vx = [xx(i6).' xx(i6).' xx(i6+1).' xx(i6+1).'];
vy = [yy(j6).' yy(j6+1).'  yy(j6+1).' yy(j6).'];
vz = [zz(k6+1).' zz(k6+1).' zz(k6+1).' zz(k6+1).'];
vx  = reshape(vx.',4*len,1);    
vy  = reshape(vy.',4*len,1);    
vz  = reshape(vz.',4*len,1);    
v6 = [vx vy vz];
f6 = [(1:4:4*(len-1)+1).' (2:4:4*(len-1)+2).' (3:4:4*(len-1)+3).' (4:4:4*(len-1)+4).'];

IND = sub2ind(siz,i6+1, j6+1, k6+1);
RGB6 = matcol(m3d(IND),:);

vertices = [v1;v2;v3;v4;v5;v6];

faces = f1; 
faces = [faces; f2+ max(max(faces))];
faces = [faces; f3+ max(max(faces))];
faces = [faces; f4+ max(max(faces))];
faces = [faces; f5+ max(max(faces))];
faces = [faces; f6+ max(max(faces))];
RGB  = [RGB1;RGB2;RGB3;RGB4;RGB5;RGB6];

if ~isempty(vertices)
    patch('Vertices',vertices,'Faces',faces,'facevertexcdata',RGB,'FaceColor','flat');
end

% display bricks with zero thickness
if exist('number_of_bricks','var')
    for ind = 1:number_of_bricks

        mtype = bricks(ind).material_type;
        
        % convert coordinates to node indices on the FDTD grid
        ni = get_node_indices(bricks(ind), fdtd_domain);
        is = ni.is; js = ni.js; ks = ni.ks; 
        ie = ni.ie; je = ni.je; ke = ni.ke; 
        
        % find the zero thickness bricks
        if (is == ie)
            len = (je-js) * (ke-ks);
            vx = zeros(len, 4); vy = zeros(len, 4); vz = zeros(len, 4);
            tind = 0;
            for mj=js:je-1
                for mk=ks:ke-1
                    tind = tind+1;
                    vx(tind,:) = [xx(is) xx(is) xx(is) xx(is)];
                    vy(tind,:) = [yy(mj) yy(mj+1) yy(mj+1) yy(mj)];
                    vz(tind,:) = [zz(mk) zz(mk) zz(mk+1) zz(mk+1)];                    
                end
            end
        end
        if (js == je)
            len = (ie-is) * (ke-ks);
            vx = zeros(len, 4); vy = zeros(len, 4); vz = zeros(len, 4);
            tind = 0;
            for mi=is:ie-1
                for mk=ks:ke-1
                    tind = tind+1;
                    vx(tind,:) = [xx(mi) xx(mi) xx(mi+1) xx(mi+1)];                    
                    vy(tind,:) = [yy(js) yy(js) yy(js) yy(js)];
                    vz(tind,:) = [zz(mk) zz(mk+1) zz(mk+1) zz(mk)];
                end
            end
        end
        if (ks == ke)
            len = (ie-is) * (je-js);
            vx = zeros(len, 4); vy = zeros(len, 4); vz = zeros(len, 4);
            tind = 0;
            for mi=is:ie-1
                for mj=js:je-1
                    tind = tind+1;
                    vx(tind,:) = [xx(mi) xx(mi+1) xx(mi+1) xx(mi)];
                    vy(tind,:) = [yy(mj) yy(mj) yy(mj+1) yy(mj+1)];                    
                    vz(tind,:) = [zz(ks) zz(ks) zz(ks) zz(ks)];
                end
            end
        end
            vx  = reshape(vx.',4*len,1);    
            vy  = reshape(vy.',4*len,1);    
            vz  = reshape(vz.',4*len,1);    
            v = [vx vy vz];
            f = [(1:4:4*(len-1)+1).' (2:4:4*(len-1)+2).' (3:4:4*(len-1)+3).' (4:4:4*(len-1)+4).'];
            RGB = matcol(mtype,:);
            patch('Vertices',v,'Faces',f,'facevertexcdata',RGB,'FaceColor','flat');
    end
end
 
display_3D_geometry;

% enable the camera toolbar
cameratoolbar(figure_number);
axis equal;
axis off;
view(20,30);

clear m3d xx yy zz vx vy vz v;
clear i1 i2 i3 i4 i5 i6 j1 j2 j3 j4 j5 j6 k1 k2 k3 k4 k5 k6;
clear v1 v2 v3 v4 v5 v6 f1 f2 f3 f4 f5 f6;
clear RGB1 RGB2 RGB3 RGB4 RGB5 RGB6 RGB;
