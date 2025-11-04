
% display spheres
for indi=1:number_of_spheres
    [X, Y, Z] = sphere(36);
    X = X * spheres(indi).radius; 
    Y = Y * spheres(indi).radius; 
    Z = Z * spheres(indi).radius; 
    X = X + spheres(indi).center_x; 
    Y = Y + spheres(indi).center_y; 
    Z = Z + spheres(indi).center_z; 
    RGB = material_types(spheres(indi).material_type).color;
    patch(surf2patch(X,Y,Z,Z),'FaceColor','none', ...
          'EdgeColor','w');
end

% display bricks
for indi=1:number_of_bricks
     vertices=[0 0 0;1 0 0 ;1 1 0;0 1 0; 0 0 1;1 0 1;1 1 1; 0 1 1 ];
     faces=[1 2 6 5;2 3 7 6;3 4 8 7;4 1 5 8;1 2 3 4;5 6 7 8 ];
     vertices(:,1)=(vertices(:,1) * (bricks(indi).max_x - bricks(indi).min_x))+bricks(indi).min_x;
     vertices(:,2)=(vertices(:,2) * (bricks(indi).max_y - bricks(indi).min_y))+bricks(indi).min_y;
     vertices(:,3)=(vertices(:,3) * (bricks(indi).max_z - bricks(indi).min_z))+bricks(indi).min_z;
     RGB = material_types(bricks(indi).material_type).color;
     patch('Vertices',vertices,'Faces',faces,...
          'FaceColor','none','edgecolor','w');
end
  

 mx = fdtd_domain.min_x;
 my = fdtd_domain.min_y;
 mz = fdtd_domain.min_z;
 
% draw axes
 line([mx 10*dx+mx],[my my], [mz mz],'color','b','linewidth',3);
 text (13*dx+mx,my,mz,'x','fontsize',12);
 line([mx mx],[my 10*dy+my], [mz mz],'color','b','linewidth',3);
 text (mx,13*dy+my,mz,'y','fontsize',12);
 line([mx mx],[my my], [mz 10*dz+mz],'color','b','linewidth',3);
 text (mx,my,13*dz+mz,'z','fontsize',12);

 
 lx = nx*dx+mx;
 ly = ny*dy+my;
 lz = nz*dz+mz;
 % draw computational space boundaries
 line([mx lx],[my my], [mz mz],'color','r','linewidth',1);
 line([mx lx],[ly ly], [mz mz],'color','r','linewidth',1);
 line([mx lx],[ly ly], [lz lz],'color','r','linewidth',1);
 line([mx lx],[my my], [lz lz],'color','r','linewidth',1);

 line([mx mx],[my ly], [mz mz],'color','r','linewidth',1);
 line([mx mx],[my ly], [lz lz],'color','r','linewidth',1);
 line([lx lx],[my ly], [mz mz],'color','r','linewidth',1);
 line([lx lx],[my ly], [lz lz],'color','r','linewidth',1);

 line([mx mx],[my my], [mz lz],'color','r','linewidth',1);
 line([mx mx],[ly ly], [mz lz],'color','r','linewidth',1);
 line([lx lx],[my my], [mz lz],'color','r','linewidth',1);
 line([lx lx],[ly ly], [mz lz],'color','r','linewidth',1);

% display thin wires
if exist('thin_wires','var')
    for indi=1:size(thin_wires,2)
        o = thin_wires(indi);
        ix = o.min_x; iy = o.min_y; iz = o.min_z; 
        ax = o.max_x; ay = o.max_y; az = o.max_z;
        line([ix ax],[iy ay],[iz az],'color',[0.8 0.3 0.1],'linewidth',2);
    end
end
