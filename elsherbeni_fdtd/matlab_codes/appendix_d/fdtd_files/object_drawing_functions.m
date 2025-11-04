function object_drawing_functions(obj)

if strcmp(obj.obj_type, 'voltage_source')
    draw_a_voltage_source(obj);
end

if strcmp(obj.obj_type, 'current_source')
    draw_a_current_source(obj);
end

if strcmp(obj.obj_type, 'resistor')
    draw_a_resistor(obj);
end

if strcmp(obj.obj_type, 'capacitor')
    draw_a_capacitor(obj);
end

if strcmp(obj.obj_type, 'inductor')
    draw_an_inductor(obj);
end

if strcmp(obj.obj_type, 'diode')
    draw_a_diode(obj);
end

if strcmp(obj.obj_type, 'sampled_voltage')
    draw_a_sampled_voltage(obj);
end

if strcmp(obj.obj_type, 'sampled_current')
    draw_a_sampled_current(obj);
end

if strcmp(obj.obj_type, 'thin_wire')
    draw_a_thin_wire(obj);
end

if strcmp(obj.obj_type, 'sampled_electric_field')
    draw_a_sampled_electric_field(obj);
end

if strcmp(obj.obj_type, 'sampled_magnetic_field')
    draw_a_sampled_magnetic_field(obj);
end

if strcmp(obj.obj_type, 'axis_outside_domain')
    draw_axis_outside_domain(obj);
end

if strcmp(obj.obj_type, 'axis_at_origin')
    draw_axis_at_origin(obj);
end

%===================================================================
function draw_a_current_source(obj)

ix = obj.min_x; iy = obj.min_y; iz = obj.min_z; 
ax = obj.max_x; ay = obj.max_y; az = obj.max_z;

switch (obj.direction)
case ('xn')
    vx = [ax ix]; vy = 0.5*[(iy+ay) (iy+ay)]; vz = 0.5*[(iz+az) (iz+az)];
    len = ax - ix;
    ir = len/20;
    [v1, f1] = a_line([ix ix+0.2*len], vy, vz, ir);
    [v2, f2] = a_line([ix+0.8*len ax], vy, vz, ir); 
    if ((az-iz)<(ay-iy)) dn = 'z'; else dn = 'y'; end
    cx = (vx(1)+vx(2))/2;
    cy = (vy(1)+vy(2))/2;
    cz = (vz(1)+vz(2))/2;
    [v3, f3] = a_circle(cx, cy, cz, len*0.3, dn, ir);    
    [v4, f4] = a_vector([cx+len/10 cx-len/10], [cy cy], [cz cz], len/16, len/8, 0.5);
case ('xp')
    vx = [ix ax]; vy = 0.5*[(iy+ay) (iy+ay)]; vz = 0.5*[(iz+az) (iz+az)];
    len = ax - ix;
    ir = len/20;
    [v1, f1] = a_line([ix ix+0.2*len], vy, vz, ir);
    [v2, f2] = a_line([ix+0.8*len ax], vy, vz, ir); 
    if ((az-iz)<(ay-iy)) dn = 'z'; else dn = 'y'; end
    cx = (vx(1)+vx(2))/2;
    cy = (vy(1)+vy(2))/2;
    cz = (vz(1)+vz(2))/2;
    [v3, f3] = a_circle(cx, cy, cz, len*0.3, dn, ir);    
    [v4, f4] = a_vector([cx-len/10 cx+len/10], [cy cy], [cz cz], len/16, len/8, 0.5);
case ('yn')
    vx = 0.5*[(ix+ax) (ix+ax)]; vy = [ay iy]; vz = 0.5*[(iz+az) (iz+az)];
    len = ay - iy;
    ir = len/20;
    [v1, f1] = a_line(vx, [iy iy+0.2*len], vz, ir);
    [v2, f2] = a_line(vx, [iy+0.8*len ay], vz, ir); 
    if ((ax-ix)<(az-iz)) dn = 'x'; else dn = 'z'; end
    cx = (vx(1)+vx(2))/2;
    cy = (vy(1)+vy(2))/2;
    cz = (vz(1)+vz(2))/2;
    [v3, f3] = a_circle(cx, cy, cz, len*0.3, dn, ir);    
    [v4, f4] = a_vector([cx cx],[cy+len/10 cy-len/10], [cz cz], len/16, len/8, 0.5);
case ('yp')
    vx = 0.5*[(ix+ax) (ix+ax)]; vy = [iy ay]; vz = 0.5*[(iz+az) (iz+az)];
    len = ay - iy;
    ir = len/20;
    [v1, f1] = a_line(vx, [iy iy+0.2*len], vz, ir);
    [v2, f2] = a_line(vx, [iy+0.8*len ay], vz, ir); 
    if ((ax-ix)<(az-iz)) dn = 'x'; else dn = 'z'; end
    cx = (vx(1)+vx(2))/2;
    cy = (vy(1)+vy(2))/2;
    cz = (vz(1)+vz(2))/2;
    [v3, f3] = a_circle(cx, cy, cz, len*0.3, dn, ir);    
    [v4, f4] = a_vector([cx cx],[cy-len/10 cy+len/10], [cz cz], len/16, len/8, 0.5);
case ('zn')
    vx = 0.5*[(ix+ax) (ix+ax)]; vy = 0.5*[(iy+ay) (iy+ay)]; vz = [az iz];
    len = az - iz;
    ir = len/20;
    [v1, f1] = a_line(vx, vy, [iz iz+0.2*len], ir);
    [v2, f2] = a_line(vx, vy, [iz+0.8*len az], ir); 
    if ((ax-ix)<(ay-iy)) dn = 'x'; else dn = 'y'; end
    cx = (vx(1)+vx(2))/2;
    cy = (vy(1)+vy(2))/2;
    cz = (vz(1)+vz(2))/2;
    [v3, f3] = a_circle(cx, cy, cz, len*0.3, dn, ir);    
    [v4, f4] = a_vector([cx cx], [cy cy], [cz+len/10 cz-len/10], len/16, len/8, 0.5);
case ('zp')
    vx = 0.5*[(ix+ax) (ix+ax)]; vy = 0.5*[(iy+ay) (iy+ay)]; vz = [iz az];
    len = az - iz;
    ir = len/20;
    [v1, f1] = a_line(vx, vy, [iz iz+0.2*len], ir);
    [v2, f2] = a_line(vx, vy, [iz+0.8*len az], ir); 
    if ((ax-ix)<(ay-iy)) dn = 'x'; else dn = 'y'; end
    cx = (vx(1)+vx(2))/2;
    cy = (vy(1)+vy(2))/2;
    cz = (vz(1)+vz(2))/2;
    [v3, f3] = a_circle(cx, cy, cz, len*0.3, dn, ir);    
    [v4, f4] = a_vector([cx cx], [cy cy], [cz-len/10 cz+len/10], len/16, len/8, 0.5);
end
v = [v1; v2; v3; v4];
f = f1;
f = [f; max(max(f)) + f2];
f = [f; max(max(f)) + f3];
f = [f; max(max(f)) + f4];
patch('Vertices',v,'Faces',f,'facecolor','m', ...
    'facealpha',1,'edgecolor','none');

% draw boundaries
v = [0 0 0; 1 0 0; 1 1 0; 0 1 0; 0 0 1; 1 0 1; 1 1 1; 0 1 1];
vx = v(:,1) * (ax - ix) + ix;
vy = v(:,2) * (ay - iy) + iy;
vz = v(:,3) * (az - iz) + iz;
f = [1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8; 1 2 3 4; 5 6 7 8];
patch('Vertices',[vx vy vz],'Faces',f,'facecolor','none','LineWidth',1, ...
    'facealpha',0.3,'edgecolor','m','linestyle','-');

% show label
text(ix,iy,iz,obj.obj_id,'VerticalAlignment','bottom', ...
    'BackgroundColor','w','edgecolor','m','fontsize',12);

%===================================================================
function draw_a_voltage_source(obj)

ix = obj.min_x; iy = obj.min_y; iz = obj.min_z; 
ax = obj.max_x; ay = obj.max_y; az = obj.max_z;

switch (obj.direction)
case ('xn')
    vx = [ax ix]; vy = 0.5*[(iy+ay) (iy+ay)]; vz = 0.5*[(iz+az) (iz+az)];
    len = ax - ix;
    ir = len/30;
    [v1, f1] = a_line([ix ix+0.2*len], vy, vz, ir);
    [v2, f2] = a_line([ix+0.8*len ax], vy, vz, ir); 
    if ((az-iz)<(ay-iy)) dn = 'z'; else dn = 'y'; end
    cx = (vx(1)+vx(2))/2;
    cy = (vy(1)+vy(2))/2;
    cz = (vz(1)+vz(2))/2;
    [v3, f3] = a_circle(cx, cy, cz, len*0.3, dn, ir);    
    [v4, f4] = a_plus(cx-len/10, cy, cz, len/5, dn, ir);
    [v5, f5] = a_minus(cx+len/10, cy, cz, len/5, dn, ir, 'x');
case ('xp')
    vx = [ix ax]; vy = 0.5*[(iy+ay) (iy+ay)]; vz = 0.5*[(iz+az) (iz+az)];
    len = ax - ix;
    ir = len/30;
    [v1, f1] = a_line([ix ix+0.2*len], vy, vz, ir);
    [v2, f2] = a_line([ix+0.8*len ax], vy, vz, ir); 
    if ((az-iz)<(ay-iy)) dn = 'z'; else dn = 'y'; end
    cx = (vx(1)+vx(2))/2;
    cy = (vy(1)+vy(2))/2;
    cz = (vz(1)+vz(2))/2;
    [v3, f3] = a_circle(cx, cy, cz, len*0.3, dn, ir);    
    [v4, f4] = a_plus(cx+len/10, cy, cz, len/5, dn, ir);
    [v5, f5] = a_minus(cx-len/10, cy, cz, len/5, dn, ir, 'x');
case ('yn')
    vx = 0.5*[(ix+ax) (ix+ax)]; vy = [ay iy]; vz = 0.5*[(iz+az) (iz+az)];
    len = ay - iy;
    ir = len/30;
    [v1, f1] = a_line(vx, [iy iy+0.2*len], vz, ir);
    [v2, f2] = a_line(vx, [iy+0.8*len ay], vz, ir); 
    if ((ax-ix)<(az-iz)) dn = 'x'; else dn = 'z'; end
    cx = (vx(1)+vx(2))/2;
    cy = (vy(1)+vy(2))/2;
    cz = (vz(1)+vz(2))/2;
    [v3, f3] = a_circle(cx, cy, cz, len*0.3, dn, ir);    
    [v4, f4] = a_plus(cx, cy-len/10, cz, len/5, dn, ir);
    [v5, f5] = a_minus(cx, cy+len/10, cz, len/5, dn, ir, 'y');
case ('yp')
    vx = 0.5*[(ix+ax) (ix+ax)]; vy = [iy ay]; vz = 0.5*[(iz+az) (iz+az)];
    len = ay - iy;
    ir = len/30;
    [v1, f1] = a_line(vx, [iy iy+0.2*len], vz, ir);
    [v2, f2] = a_line(vx, [iy+0.8*len ay], vz, ir); 
    if ((ax-ix)<(az-iz)) dn = 'x'; else dn = 'z'; end
    cx = (vx(1)+vx(2))/2;
    cy = (vy(1)+vy(2))/2;
    cz = (vz(1)+vz(2))/2;
    [v3, f3] = a_circle(cx, cy, cz, len*0.3, dn, ir);    
    [v4, f4] = a_plus(cx, cy+len/10, cz, len/5, dn, ir);
    [v5, f5] = a_minus(cx, cy-len/10, cz, len/5, dn, ir, 'y');    
case ('zn')
    vx = 0.5*[(ix+ax) (ix+ax)]; vy = 0.5*[(iy+ay) (iy+ay)]; vz = [az iz];
    len = az - iz;
    ir = len/30;
    [v1, f1] = a_line(vx, vy, [iz iz+0.2*len], ir);
    [v2, f2] = a_line(vx, vy, [iz+0.8*len az], ir); 
    if ((ax-ix)<(ay-iy)) dn = 'x'; else dn = 'y'; end
    cx = (vx(1)+vx(2))/2;
    cy = (vy(1)+vy(2))/2;
    cz = (vz(1)+vz(2))/2;
    [v3, f3] = a_circle(cx, cy, cz, len*0.3, dn, ir);    
    [v4, f4] = a_plus(cx, cy, cz-len/10, len/5, dn, ir);
    [v5, f5] = a_minus(cx, cy, cz+len/10, len/5, dn, ir, 'z');
case ('zp')
    vx = 0.5*[(ix+ax) (ix+ax)]; vy = 0.5*[(iy+ay) (iy+ay)]; vz = [iz az];
    len = az - iz;
    ir = len/30;
    [v1, f1] = a_line(vx, vy, [iz iz+0.2*len], ir);
    [v2, f2] = a_line(vx, vy, [iz+0.8*len az], ir); 
    if ((ax-ix)<(ay-iy)) dn = 'x'; else dn = 'y'; end
    cx = (vx(1)+vx(2))/2;
    cy = (vy(1)+vy(2))/2;
    cz = (vz(1)+vz(2))/2;
    [v3, f3] = a_circle(cx, cy, cz, len*0.3, dn, ir);    
    [v4, f4] = a_plus(cx, cy, cz+len/10, len/5, dn, ir);
    [v5, f5] = a_minus(cx, cy, cz-len/10, len/5, dn, ir, 'z');
end
v = [v1; v2; v3; v4; v5];
f = f1;
f = [f; max(max(f)) + f2];
f = [f; max(max(f)) + f3];
f = [f; max(max(f)) + f4];
f = [f; max(max(f)) + f5];
patch('Vertices',v,'Faces',f,'facecolor','g', ...
    'facealpha',1,'edgecolor','none');

% draw boundaries
v = [0 0 0; 1 0 0; 1 1 0; 0 1 0; 0 0 1; 1 0 1; 1 1 1; 0 1 1];
vx = v(:,1) * (ax - ix) + ix;
vy = v(:,2) * (ay - iy) + iy;
vz = v(:,3) * (az - iz) + iz;
f = [1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8; 1 2 3 4; 5 6 7 8];
patch('Vertices',[vx vy vz],'Faces',f,'facecolor','none','LineWidth',1, ...
    'facealpha',0.3,'edgecolor','g','linestyle','-');

% show label
text(ix,iy,iz,obj.obj_id,'VerticalAlignment','bottom', ...
    'BackgroundColor','w','edgecolor','g','fontsize',12);

%===================================================================
function draw_a_resistor(obj)

ix = obj.min_x; iy = obj.min_y; iz = obj.min_z; 
ax = obj.max_x; ay = obj.max_y; az = obj.max_z;

t = 0.1 * ones(1,11);
[X,Y,Z] = cylinder(t,18);
X(5,:) = X(5,:) + 0.2;
X(7,:) = X(7,:) - 0.2;

switch (obj.direction)
case ('x')
    len = ax - ix;
    vx = ix; vy = 0.5*(iy+ay); vz = 0.5*(iz+az);
    X = X*len; Y = Y*len; Z = Z*len;    
    if ((az-iz)<(ay-iy)) 
        [f,v] = surf2patch(Z+vx,X+vy,Y+vz,Z);       
    else
        [f,v] = surf2patch(Z+vx,Y+vy,X+vz,Z);       
    end
case ('y')
    len = ay - iy;
    vx = 0.5*(ix+ax); vy = iy; vz = 0.5*(iz+az);
    X = X*len; Y = Y*len; Z = Z*len;    
    if ((ax-ix)<(az-iz))
        [f,v] = surf2patch(Y+vx,Z+vy,X+vz,Z);       
    else
        [f,v] = surf2patch(X+vx,Z+vy,Y+vz,Z);       
    end
case ('z')
    len = az - iz;
    vy = 0.5*(iy+ay); vx = 0.5*(ix+ax); vz = iz;
    X = X*len; Y = Y*len; Z = Z*len;    
    if ((ax-ix)<(ay-iy)) 
        [f,v] = surf2patch(Y+vx,X+vy,Z+vz,Z);       
    else
        [f,v] = surf2patch(X+vx,Y+vy,Z+vz,Z);       
    end
end
patch('Vertices',v,'Faces',f,'facecolor','k', ...
    'facealpha',1,'edgecolor','none');

% draw boundaries
v = [0 0 0; 1 0 0; 1 1 0; 0 1 0; 0 0 1; 1 0 1; 1 1 1; 0 1 1];
vx = v(:,1) * (ax - ix) + ix;
vy = v(:,2) * (ay - iy) + iy;
vz = v(:,3) * (az - iz) + iz;
f = [1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8; 1 2 3 4; 5 6 7 8];
patch('Vertices',[vx vy vz],'Faces',f,'facecolor','none','LineWidth',1, ...
    'facealpha',0.3,'edgecolor','k','linestyle','-');

% show label
text(ix,iy,iz,obj.obj_id,'VerticalAlignment','bottom', ...
    'BackgroundColor','w','edgecolor','k','fontsize',12);

%===================================================================
function draw_a_capacitor(obj)

ix = obj.min_x; iy = obj.min_y; iz = obj.min_z; 
ax = obj.max_x; ay = obj.max_y; az = obj.max_z;

t = 0.1 * ones(1,11);
[X,Y,Z] = cylinder(t,18);
X(5,:) = X(5,:) + 0.2;
X(7,:) = X(7,:) - 0.2;

switch (obj.direction)
case ('x')
    len = ax - ix;
    vx = ix; vy = 0.5*(iy+ay); vz = 0.5*(iz+az);
    ir = len/20;
    [v1, f1] = a_line([ix ix+0.4*len], [vy vy], [vz vz], ir);
    [v2, f2] = a_line([ix+0.6*len ax], [vy vy], [vz vz], ir);
    [v3, f3] = a_line([ix+0.4*len ix+0.42*len], [vy vy], [vz vz], 4*ir);
    [v4, f4] = a_line([ix+0.58*len ix+0.6*len], [vy vy], [vz vz], 4*ir);
case ('y')
    len = ay - iy;
    vx = 0.5*(ix+ax); vy = iy; vz = 0.5*(iz+az);
    ir = len/20;
    [v1, f1] = a_line([vx vx], [iy iy+0.4*len], [vz vz], ir);
    [v2, f2] = a_line([vx vx], [iy+0.6*len ay], [vz vz], ir);
    [v3, f3] = a_line([vx vx], [iy+0.4*len iy+0.42*len], [vz vz], 4*ir);
    [v4, f4] = a_line([vx vx], [iy+0.58*len iy+0.6*len], [vz vz], 4*ir);
case ('z')
    len = az - iz;
    vy = 0.5*(iy+ay); vx = 0.5*(ix+ax); vz = iz;
    ir = len/20;
    [v1, f1] = a_line([vx vx], [vy vy], [iz iz+0.4*len], ir);
    [v2, f2] = a_line([vx vx], [vy vy], [iz+0.6*len az], ir);
    [v3, f3] = a_line([vx vx], [vy vy], [iz+0.4*len iz+0.42*len], 4*ir);
    [v4, f4] = a_line([vx vx], [vy vy], [iz+0.58*len iz+0.6*len], 4*ir);
end
v = [v1; v2; v3; v4];
f = f1;
f = [f; max(max(f)) + f2];
f = [f; max(max(f)) + f3];
f = [f; max(max(f)) + f4];
patch('Vertices',v,'Faces',f,'facecolor','r', ...
    'facealpha',1,'edgecolor','none');

% draw boundaries
v = [0 0 0; 1 0 0; 1 1 0; 0 1 0; 0 0 1; 1 0 1; 1 1 1; 0 1 1];
vx = v(:,1) * (ax - ix) + ix;
vy = v(:,2) * (ay - iy) + iy;
vz = v(:,3) * (az - iz) + iz;
f = [1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8; 1 2 3 4; 5 6 7 8];
patch('Vertices',[vx vy vz],'Faces',f,'facecolor','none','LineWidth',1, ...
    'facealpha',0.3,'edgecolor','r','linestyle','-');

% show label
text(ix,iy,iz,obj.obj_id,'VerticalAlignment','bottom', ...
    'BackgroundColor','w','edgecolor','r','fontsize',12);

%===================================================================
function draw_an_inductor(obj)

ix = obj.min_x; iy = obj.min_y; iz = obj.min_z; 
ax = obj.max_x; ay = obj.max_y; az = obj.max_z;

t = 0.05 * ones(1,41);
[X,Y,Z] = cylinder(t,18);
t = 0.2*abs(sin(2*pi*[0:22]/22));
for mi = 1:23 
    X(mi+8,:) = X(mi+8,:) + t(mi);
end

switch (obj.direction)
case ('x')
    len = ax - ix;
    vx = ix; vy = 0.5*(iy+ay); vz = 0.5*(iz+az);
    X = X*len; Y = Y*len; Z = Z*len;    
    if ((az-iz)<(ay-iy)) 
        [f,v] = surf2patch(Z+vx,X+vy,Y+vz,Z);       
    else
        [f,v] = surf2patch(Z+vx,Y+vy,X+vz,Z);       
    end
case ('y')
    len = ay - iy;
    vx = 0.5*(ix+ax); vy = iy; vz = 0.5*(iz+az);
    X = X*len; Y = Y*len; Z = Z*len;    
    if ((ax-ix)<(az-iz))
        [f,v] = surf2patch(Y+vx,Z+vy,X+vz,Z);       
    else
        [f,v] = surf2patch(X+vx,Z+vy,Y+vz,Z);       
    end
case ('z')
    len = az - iz;
    vy = 0.5*(iy+ay); vx = 0.5*(ix+ax); vz = iz;
    X = X*len; Y = Y*len; Z = Z*len;    
    if ((ax-ix)<(ay-iy)) 
        [f,v] = surf2patch(Y+vx,X+vy,Z+vz,Z);       
    else
        [f,v] = surf2patch(X+vx,Y+vy,Z+vz,Z);       
    end
end
patch('Vertices',v,'Faces',f,'facecolor','b', ...
    'facealpha',1,'edgecolor','none');

% draw boundaries
v = [0 0 0; 1 0 0; 1 1 0; 0 1 0; 0 0 1; 1 0 1; 1 1 1; 0 1 1];
vx = v(:,1) * (ax - ix) + ix;
vy = v(:,2) * (ay - iy) + iy;
vz = v(:,3) * (az - iz) + iz;
f = [1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8; 1 2 3 4; 5 6 7 8];
patch('Vertices',[vx vy vz],'Faces',f,'facecolor','none','LineWidth',1, ...
    'facealpha',0.3,'edgecolor','b','linestyle','-');

% show label
text(ix,iy,iz,obj.obj_id,'VerticalAlignment','bottom', ...
    'BackgroundColor','w','edgecolor','b','fontsize',12);

%===================================================================
function draw_a_thin_wire(obj)

ix = obj.min_x; iy = obj.min_y; iz = obj.min_z; 
ax = obj.max_x; ay = obj.max_y; az = obj.max_z;

ir = obj.radius;
switch (obj.direction)
    case 'x'
        [v, f] = a_line([ix ax],[iy iy],[iz iz],ir);
    case 'y'
        [v, f] = a_line([ix ix],[iy ay],[iz iz],ir);
    case 'z'
        [v, f] = a_line([ix ix],[iy iy],[iz az],ir);
end
patch('Vertices',v,'Faces',f,'facecolor','r', ...
    'facealpha',1,'edgecolor','none');

% show label
text(ix,iy,iz,obj.obj_id,'VerticalAlignment','bottom', ...
    'BackgroundColor','w','edgecolor','r','fontsize',12);

%===================================================================
function draw_a_sampled_electric_field(obj)

ix = obj.x; iy = obj.y; iz = obj.z; 
dx = obj.dx; dy = obj.dy; dz = obj.dz;

ir = max([dx dy dz])/5;

[v1, f1] = a_line([ix ix],[iy iy],[iz-2*ir iz+2*ir],ir/2);
[v2, f2] = a_line([ix ix+2*ir],[iy iy],[iz-2*ir iz-2*ir],ir/2);
[v3, f3] = a_line([ix ix+2*ir],[iy iy],[iz iz],ir/2);
[v4, f4] = a_line([ix ix+2*ir],[iy iy],[iz+2*ir iz+2*ir],ir/2);

v = [v1; v2; v3; v4];
f = f1;
f = [f; max(max(f)) + f2];
f = [f; max(max(f)) + f3];
f = [f; max(max(f)) + f4];

patch('Vertices',v,'Faces',f,'facecolor','r', ...
    'facealpha',1,'edgecolor','none');

% show label
text(ix,iy,iz,obj.obj_id,'VerticalAlignment','bottom', ...
    'BackgroundColor','w','edgecolor','r','fontsize',12);

%===================================================================
function draw_a_sampled_magnetic_field(obj)

ix = obj.x; iy = obj.y; iz = obj.z; 
dx = obj.dx; dy = obj.dy; dz = obj.dz;

ir = max([dx dy dz])/5;

[v1, f1] = a_line([ix ix],[iy iy],[iz-2*ir iz+2*ir],ir/2);
[v2, f2] = a_line([ix+2*ir ix+2*ir],[iy iy],[iz-2*ir iz+2*ir],ir/2);
[v3, f3] = a_line([ix ix+2*ir],[iy iy],[iz iz],ir/2);

v = [v1; v2; v3];
f = f1;
f = [f; max(max(f)) + f2];
f = [f; max(max(f)) + f3];

patch('Vertices',v,'Faces',f,'facecolor','b', ...
    'facealpha',1,'edgecolor','none');

% show label
text(ix,iy,iz,obj.obj_id,'VerticalAlignment','bottom', ...
    'BackgroundColor','w','edgecolor','b','fontsize',12);

%===================================================================
function draw_a_diode(obj)

ix = obj.min_x; iy = obj.min_y; iz = obj.min_z; 
ax = obj.max_x; ay = obj.max_y; az = obj.max_z;

t = 0.1 * ones(1,11);
[X,Y,Z] = cylinder(t,18);
X(5,:) = X(5,:) + 0.2;
X(7,:) = X(7,:) - 0.2;

switch (obj.direction)
case ('xn')
    len = ax - ix;
    vx = ix; vy = 0.5*(iy+ay); vz = 0.5*(iz+az);
    ir = len/20;
    [v1, f1] = a_line([ix ax], [vy vy], [vz vz], ir);
    [v2, f2] = a_vector([ax ix+len*0.25], [vy vy], [vz vz], ir, 4*ir, 0.5);
    [v3, f3] = a_line([ix+0.34*len ix+0.36*len], [vy vy], [vz vz], 4*ir);
case ('xp')
    len = ax - ix;
    vx = ix; vy = 0.5*(iy+ay); vz = 0.5*(iz+az);
    ir = len/20;
    [v1, f1] = a_line([ix ax], [vy vy], [vz vz], ir);
    [v2, f2] = a_vector([ix ax-len*0.25], [vy vy], [vz vz], ir, 4*ir, 0.5);
    [v3, f3] = a_line([ix+0.64*len ix+0.66*len], [vy vy], [vz vz], 4*ir);
case ('yn')
    len = ay - iy;
    vx = 0.5*(ix+ax); vy = iy; vz = 0.5*(iz+az);
    ir = len/20;
    [v1, f1] = a_line([vx vx], [iy ay], [vz vz], ir);
    [v2, f2] = a_vector([vx vx], [ay iy+len*0.25], [vz vz], ir, 4*ir, 0.5);
    [v3, f3] = a_line([vx vx], [iy+0.34*len iy+0.36*len], [vz vz], 4*ir);
case ('yp')
    len = ay - iy;
    vx = 0.5*(ix+ax); vy = iy; vz = 0.5*(iz+az);
    ir = len/20;
    [v1, f1] = a_line([vx vx], [iy ay], [vz vz], ir);
    [v2, f2] = a_vector([vx vx], [iy ay-len*0.25], [vz vz], ir, 4*ir, 0.5);
    [v3, f3] = a_line([vx vx], [iy+0.64*len iy+0.66*len], [vz vz], 4*ir);
case ('zn')
    len = az - iz;
    vy = 0.5*(iy+ay); vx = 0.5*(ix+ax); vz = iz;
    ir = len/20;
    [v1, f1] = a_line([vx vx], [vy vy], [iz az], ir);
    [v2, f2] = a_vector([vx vx], [vy vy], [az iz+len*0.25], ir, 4*ir, 0.5);
    [v3, f3] = a_line([vx vx], [vy vy], [iz+0.34*len iz+0.36*len], 4*ir);
case ('zp')
    len = az - iz;
    vy = 0.5*(iy+ay); vx = 0.5*(ix+ax); vz = iz;
    ir = len/20;
    [v1, f1] = a_line([vx vx], [vy vy], [iz az], ir);
    [v2, f2] = a_vector([vx vx], [vy vy], [iz az-len*0.25], ir, 4*ir, 0.5);
    [v3, f3] = a_line([vx vx], [vy vy], [iz+0.64*len iz+0.66*len], 4*ir);
end
v = [v1; v2; v3];
f = f1;
f = [f; max(max(f)) + f2];
f = [f; max(max(f)) + f3];
col = [0.5 0.5 0.5];
patch('Vertices',v,'Faces',f,'facecolor',col, ...
    'facealpha',1,'edgecolor','none');

% draw boundaries
v = [0 0 0; 1 0 0; 1 1 0; 0 1 0; 0 0 1; 1 0 1; 1 1 1; 0 1 1];
vx = v(:,1) * (ax - ix) + ix;
vy = v(:,2) * (ay - iy) + iy;
vz = v(:,3) * (az - iz) + iz;
f = [1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8; 1 2 3 4; 5 6 7 8];
patch('Vertices',[vx vy vz],'Faces',f,'facecolor','none','LineWidth',1, ...
    'facealpha',0.3,'edgecolor',col,'linestyle','-');

% show label
text(ix,iy,iz,obj.obj_id,'VerticalAlignment','bottom', ...
    'BackgroundColor','w','edgecolor',col,'fontsize',12);

%===================================================================
function draw_a_sampled_voltage(obj)

ix = obj.min_x; iy = obj.min_y; iz = obj.min_z; 
ax = obj.max_x; ay = obj.max_y; az = obj.max_z;

switch (obj.direction)
case ('xn')
    len = ax - ix;
    ir = len/20;
    if ((az-iz)<(ay-iy)) dn = 'z'; else dn = 'y'; end
    cx = (ix+ax)/2; cy = (iy+ay)/2; cz = (iz+az)/2;
    [v1, f1] = a_plus(ix, iy, iz, len/5, dn, ir);
    [v2, f2] = a_plus(ix, iy, az, len/5, dn, ir);
    [v3, f3] = a_plus(ix, ay, iz, len/5, dn, ir);
    [v4, f4] = a_plus(ix, ay, az, len/5, dn, ir);
    [v6, f6] = a_minus(ax, iy, iz, len/5, dn, ir, 'x');
    [v7, f7] = a_minus(ax, iy, az, len/5, dn, ir, 'x');
    [v8, f8] = a_minus(ax, ay, iz, len/5, dn, ir, 'x');
    [v9, f9] = a_minus(ax, ay, az, len/5, dn, ir, 'x');
case ('xp')
    len = ax - ix;
    ir = len/20;
    if ((az-iz)<(ay-iy)) dn = 'z'; else dn = 'y'; end
    cx = (ix+ax)/2; cy = (iy+ay)/2; cz = (iz+az)/2;
    [v1, f1] = a_plus(ax, iy, iz, len/5, dn, ir);
    [v2, f2] = a_plus(ax, iy, az, len/5, dn, ir);
    [v3, f3] = a_plus(ax, ay, iz, len/5, dn, ir);
    [v4, f4] = a_plus(ax, ay, az, len/5, dn, ir);
    [v6, f6] = a_minus(ix, iy, iz, len/5, dn, ir, 'x');
    [v7, f7] = a_minus(ix, iy, az, len/5, dn, ir, 'x');
    [v8, f8] = a_minus(ix, ay, iz, len/5, dn, ir, 'x');
    [v9, f9] = a_minus(ix, ay, az, len/5, dn, ir, 'x');
case ('yn')
    len = ay - iy;
    ir = len/20;
    if ((ax-ix)<(az-iz)) dn = 'x'; else dn = 'z'; end
    cx = (ix+ax)/2; cy = (iy+ay)/2; cz = (iz+az)/2;
    [v1, f1] = a_plus(ix, iy, iz, len/5, dn, ir);
    [v2, f2] = a_plus(ix, iy, az, len/5, dn, ir);
    [v3, f3] = a_plus(ax, iy, iz, len/5, dn, ir);
    [v4, f4] = a_plus(ax, iy, az, len/5, dn, ir);
    [v6, f6] = a_minus(ix, ay, iz, len/5, dn, ir, 'y');
    [v7, f7] = a_minus(ix, ay, az, len/5, dn, ir, 'y');
    [v8, f8] = a_minus(ax, ay, iz, len/5, dn, ir, 'y');
    [v9, f9] = a_minus(ax, ay, az, len/5, dn, ir, 'y');
case ('yp')
    len = ay - iy;
    ir = len/20;
    if ((ax-ix)<(az-iz)) dn = 'x'; else dn = 'z'; end
    cx = (ix+ax)/2; cy = (iy+ay)/2; cz = (iz+az)/2;
    [v1, f1] = a_plus(ix, ay, iz, len/5, dn, ir);
    [v2, f2] = a_plus(ix, ay, az, len/5, dn, ir);
    [v3, f3] = a_plus(ax, ay, iz, len/5, dn, ir);
    [v4, f4] = a_plus(ax, ay, az, len/5, dn, ir);
    [v6, f6] = a_minus(ix, iy, iz, len/5, dn, ir, 'y');
    [v7, f7] = a_minus(ix, iy, az, len/5, dn, ir, 'y');
    [v8, f8] = a_minus(ax, iy, iz, len/5, dn, ir, 'y');
    [v9, f9] = a_minus(ax, iy, az, len/5, dn, ir, 'y');
case ('zn')
    len = az - iz;
    ir = len/20;
    if ((ax-ix)<(ay-iy)) dn = 'x'; else dn = 'y'; end
    cx = (ix+ax)/2; cy = (iy+ay)/2; cz = (iz+az)/2;
    [v1, f1] = a_plus(ix, iy, iz, len/5, dn, ir);
    [v2, f2] = a_plus(ix, ay, iz, len/5, dn, ir);
    [v3, f3] = a_plus(ax, iy, iz, len/5, dn, ir);
    [v4, f4] = a_plus(ax, ay, iz, len/5, dn, ir);
    [v6, f6] = a_minus(ix, iy, az, len/5, dn, ir, 'z');
    [v7, f7] = a_minus(ix, ay, az, len/5, dn, ir, 'z');
    [v8, f8] = a_minus(ax, iy, az, len/5, dn, ir, 'z');
    [v9, f9] = a_minus(ax, ay, az, len/5, dn, ir, 'z');
case ('zp')
    len = az - iz;
    ir = len/20;
    if ((ax-ix)<(ay-iy)) dn = 'x'; else dn = 'y'; end
    cx = (ix+ax)/2; cy = (iy+ay)/2; cz = (iz+az)/2;
    [v1, f1] = a_plus(ix, iy, az, len/5, dn, ir);
    [v2, f2] = a_plus(ix, ay, az, len/5, dn, ir);
    [v3, f3] = a_plus(ax, iy, az, len/5, dn, ir);
    [v4, f4] = a_plus(ax, ay, az, len/5, dn, ir);
    [v6, f6] = a_minus(ix, iy, iz, len/5, dn, ir, 'z');
    [v7, f7] = a_minus(ix, ay, iz, len/5, dn, ir, 'z');
    [v8, f8] = a_minus(ax, iy, iz, len/5, dn, ir, 'z');
    [v9, f9] = a_minus(ax, ay, iz, len/5, dn, ir, 'z');
end
v = [v1; v2; v3; v4; v6; v7; v8; v9];
f = f1;
f = [f; max(max(f)) + f2];
f = [f; max(max(f)) + f3];
f = [f; max(max(f)) + f4];
f = [f; max(max(f)) + f6];
f = [f; max(max(f)) + f7];
f = [f; max(max(f)) + f8];
f = [f; max(max(f)) + f9];
patch('Vertices',v,'Faces',f,'facecolor','c', ...
    'facealpha',1,'edgecolor','none');

% draw boundaries
v = [0 0 0; 1 0 0; 1 1 0; 0 1 0; 0 0 1; 1 0 1; 1 1 1; 0 1 1];
vx = v(:,1) * (ax - ix) + ix;
vy = v(:,2) * (ay - iy) + iy;
vz = v(:,3) * (az - iz) + iz;
f = [1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8; 1 2 3 4; 5 6 7 8];
patch('Vertices',[vx vy vz],'Faces',f,'facecolor','none','LineWidth',1, ...
    'facealpha',0.3,'edgecolor','c','linestyle','-');

% show label
text(ix,iy,iz,obj.obj_id,'VerticalAlignment','bottom', ...
    'BackgroundColor','w','edgecolor','c','fontsize',12);

%===================================================================
function draw_a_sampled_current(obj)

ix = obj.min_x; iy = obj.min_y; iz = obj.min_z; 
ax = obj.max_x; ay = obj.max_y; az = obj.max_z;
len = max([ax-ix ay-iy az-iz]);
cx = (ix+ax)/2; cy = (iy+ay)/2; cz = (iz+az)/2;
ir = len/20;

switch (obj.direction)
case ('xn')
    [v1, f1] = a_vector([cx+3*ir cx-3*ir], [iy iy], [iz iz], ir, 2*ir, 0.5);
    [v2, f2] = a_vector([cx+3*ir cx-3*ir], [iy iy], [az az], ir, 2*ir, 0.5);
    [v3, f3] = a_vector([cx+3*ir cx-3*ir], [ay ay], [iz iz], ir, 2*ir, 0.5);
    [v4, f4] = a_vector([cx+3*ir cx-3*ir], [ay ay], [az az], ir, 2*ir, 0.5);
    [v5, f5] = a_vector([cx+3*ir cx-3*ir], [cy cy], [cz cz], ir, 2*ir, 0.5);
case ('xp')
    [v1, f1] = a_vector([cx-3*ir cx+3*ir], [iy iy], [iz iz], ir, 2*ir, 0.5);
    [v2, f2] = a_vector([cx-3*ir cx+3*ir], [iy iy], [az az], ir, 2*ir, 0.5);
    [v3, f3] = a_vector([cx-3*ir cx+3*ir], [ay ay], [iz iz], ir, 2*ir, 0.5);
    [v4, f4] = a_vector([cx-3*ir cx+3*ir], [ay ay], [az az], ir, 2*ir, 0.5);
    [v5, f5] = a_vector([cx-3*ir cx+3*ir], [cy cy], [cz cz], ir, 2*ir, 0.5);
case ('yn')
    [v1, f1] = a_vector([ix ix], [cy+3*ir cy-3*ir], [iz iz],  ir, 2*ir, 0.5);
    [v2, f2] = a_vector([ax ax], [cy+3*ir cy-3*ir], [iz iz],  ir, 2*ir, 0.5);
    [v3, f3] = a_vector([ix ix], [cy+3*ir cy-3*ir], [az az],  ir, 2*ir, 0.5);
    [v4, f4] = a_vector([ax ax], [cy+3*ir cy-3*ir], [az az],  ir, 2*ir, 0.5);
    [v5, f5] = a_vector([cx cx], [cy+3*ir cy-3*ir], [cz cz],  ir, 2*ir, 0.5);
case ('yp')
    [v1, f1] = a_vector([ix ix], [cy-3*ir cy+3*ir], [iz iz],  ir, 2*ir, 0.5);
    [v2, f2] = a_vector([ax ax], [cy-3*ir cy+3*ir], [iz iz],  ir, 2*ir, 0.5);
    [v3, f3] = a_vector([ix ix], [cy-3*ir cy+3*ir], [az az],  ir, 2*ir, 0.5);
    [v4, f4] = a_vector([ax ax], [cy-3*ir cy+3*ir], [az az],  ir, 2*ir, 0.5);
    [v5, f5] = a_vector([cx cx], [cy-3*ir cy+3*ir], [cz cz],  ir, 2*ir, 0.5);
case ('zn')
    [v1, f1] = a_vector([ix ix], [iy iy], [cz+3*ir cz-3*ir],   ir, 2*ir, 0.5);
    [v2, f2] = a_vector([ix ix], [ay ay], [cz+3*ir cz-3*ir],   ir, 2*ir, 0.5);
    [v3, f3] = a_vector([ax ax], [iy iy], [cz+3*ir cz-3*ir],   ir, 2*ir, 0.5);
    [v4, f4] = a_vector([ax ax], [ay ay], [cz+3*ir cz-3*ir],   ir, 2*ir, 0.5);
    [v5, f5] = a_vector([cx cx], [cy cy], [cz+3*ir cz-3*ir],   ir, 2*ir, 0.5);
case ('zp')
    [v1, f1] = a_vector([ix ix], [iy iy], [cz-3*ir cz+3*ir],   ir, 2*ir, 0.5);
    [v2, f2] = a_vector([ix ix], [ay ay], [cz-3*ir cz+3*ir],   ir, 2*ir, 0.5);
    [v3, f3] = a_vector([ax ax], [iy iy], [cz-3*ir cz+3*ir],   ir, 2*ir, 0.5);
    [v4, f4] = a_vector([ax ax], [ay ay], [cz-3*ir cz+3*ir],   ir, 2*ir, 0.5);
    [v5, f5] = a_vector([cx cx], [cy cy], [cz-3*ir cz+3*ir],   ir, 2*ir, 0.5);
end
v = [v1; v2; v3; v4; v5];
f = f1;
f = [f; max(max(f)) + f2];
f = [f; max(max(f)) + f3];
f = [f; max(max(f)) + f4];
f = [f; max(max(f)) + f5];
col = [1 0.6 0.2];
patch('Vertices',v,'Faces',f,'facecolor',col, ...
    'facealpha',1,'edgecolor','none');

% draw boundaries
v = [0 0 0; 1 0 0; 1 1 0; 0 1 0; 0 0 1; 1 0 1; 1 1 1; 0 1 1];
vx = v(:,1) * (ax - ix) + ix;
vy = v(:,2) * (ay - iy) + iy;
vz = v(:,3) * (az - iz) + iz;
f = [1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8; 1 2 3 4; 5 6 7 8];
patch('Vertices',[vx vy vz],'Faces',f,'facecolor','none','LineWidth',1, ...
    'facealpha',0.3,'edgecolor',col,'linestyle','-');

% show label
text(ix,iy,iz,obj.obj_id,'VerticalAlignment','bottom', ...
    'BackgroundColor','w','edgecolor',col,'fontsize',12);

%========================================================================  
function draw_axis_outside_domain(obj)

dx = obj.dx;
dy = obj.dy;
dz = obj.dz;
lx = obj.lx;
ly = obj.ly;
lz = obj.lz;
d = max([dx dy dz]);

refx = lx - 10*d;
refy = ly - 10*d;
refz = lz;

[vt1, f1] = a_vector([refx refx], [refy refy], [refz refz+5*d], 0.5*d, d, 0.3);
[vt2, f2] = a_vector([refx refx], [refy refy+5*d], [refz refz], 0.5*d, d, 0.3);
[vt3, f3] = a_vector([refx refx+5*d], [refy refy], [refz refz], 0.5*d, d, 0.3);
vertices = [vt1; vt2; vt3];

f2 = f2 + max(max(f1));
f3 = f3 + max(max(f2));
faces  = [f1;f2;f3];

col = 'g';
patch('Vertices',vertices,'Faces',faces,'facecolor',col, ...
    'facealpha',1,'edgecolor','none');

text(refx+6*d,refy,refz,'x','fontsize',12);
text(refx,refy+6*d,refz,'y','fontsize',12);
text(refx,refy,refz+6*d,'z','fontsize',12);

%========================================================================  
function draw_axis_at_origin(obj)

dx = obj.dx;
dy = obj.dy;
dz = obj.dz;
d = max([dx dy dz]);
[vt1, f1] = a_vector([0 0], [0 0], [0 5*d], 0.5*d, d, 0.3);
[vt2, f2] = a_vector([0 0], [0 5*d], [0 0], 0.5*d, d, 0.3);
[vt3, f3] = a_vector([0 5*d], [0 0], [0 0], 0.5*d, d, 0.3);
vertices = [vt1; vt2; vt3];

f2 = f2 + max(max(f1));
f3 = f3 + max(max(f2));
faces  = [f1;f2;f3];

col = 'g';
patch('Vertices',vertices,'Faces',faces,'facecolor',col, ...
    'facealpha',1,'edgecolor','none');

text(6*d,0,0,'x','fontsize',12);
text(0,6*d,0,'y','fontsize',12);
text(0,0,6*d,'z','fontsize',12);

% ===============================================================
function [vertices, faces] = a_circle(cx, cy, cz, ir, dn, lr)

% cx: x coordinates of the circle center
% cy: y coordinates of the circle center
% cz: z coordinates of the circle center
% ir: radius of the circle
% dn: normal direction
% lr: line radius

t = -lr:lr/20:lr;
[X,Y,Z] = cylinder(ir-lr+sqrt(lr^2-t.^2),72);
Z = Z*2*lr-lr;
[f1,v1] = surf2patch(X,Y,Z,Z);
[X,Y,Z] = cylinder(ir-lr-sqrt(lr^2-t.^2),72);
Z = Z*2*lr-lr;
[f2,v2] = surf2patch(X,Y,Z,Z);
faces = [f1; f2+max(max(f1))];
v = [v1;v2];

switch dn
    case 'x'
        v(:,1) = v(:,1) + cy;
        v(:,2) = v(:,2) + cz;
        v(:,3) = v(:,3) + cx;
        vertices = [v(:,3) v(:,1) v(:,2)];
    case 'y'
        v(:,1) = v(:,1) + cz;
        v(:,2) = v(:,2) + cx;
        v(:,3) = v(:,3) + cy;
        vertices = [v(:,2) v(:,3) v(:,1)];
    case 'z'
        v(:,1) = v(:,1) + cx;
        v(:,2) = v(:,2) + cy;
        v(:,3) = v(:,3) + cz;
        vertices = [v(:,1) v(:,2) v(:,3)];
end

% ===============================================================
function [vertices, faces] = a_plus(cx, cy, cz, len, dn, lr)

% cx: x coordinates of the plus center
% cy: y coordinates of the plus center
% cz: z coordinates of the plus center
% len: size of the plus
% dn: normal direction
% lr: line radius

switch dn
    case 'x'
        [v1, f1] = a_line([cx cx], [cy-len/2 cy+len/2], [cz cz], lr);
        [v2, f2] = a_line([cx cx], [cy cy], [cz-len/2 cz+len/2], lr);
    case 'y'
        [v1, f1] = a_line([cx-len/2 cx+len/2], [cy cy], [cz cz], lr);
        [v2, f2] = a_line([cx cx], [cy cy], [cz-len/2 cz+len/2], lr);
    case 'z'
        [v1, f1] = a_line([cx cx], [cy-len/2 cy+len/2], [cz cz], lr);
        [v2, f2] = a_line([cx-len/2 cx+len/2], [cy cy], [cz cz], lr);
end
vertices = [v1; v2];
faces = [f1; max(max(f1)) + f2]; 

%========================================================================  
function [vertices, faces] = a_minus(cx, cy, cz, len, dn, lr, dir)

% cx: x coordinates of the plus center
% cy: y coordinates of the plus center
% cz: z coordinates of the plus center
% len: size of the plus
% dn: normal direction
% lr: line radius
% dir: text direction

switch dn
    case 'x'
        if dir=='z'
            [v, f] = a_line([cx cx], [cy-len/2 cy+len/2], [cz cz], lr);
        else
            [v, f] = a_line([cx cx], [cy cy], [cz-len/2 cz+len/2], lr);
        end
    case 'y'
        if dir=='z'
            [v, f] = a_line([cx-len/2 cx+len/2], [cy cy], [cz cz], lr);
        else
            [v, f] = a_line([cx cx], [cy cy], [cz-len/2 cz+len/2], lr);
        end
    case 'z'
        if dir=='x'
            [v, f] = a_line([cx cx], [cy-len/2 cy+len/2], [cz cz], lr);
        else
            [v, f] = a_line([cx-len/2 cx+len/2], [cy cy], [cz cz], lr);
        end
end
vertices = v;
faces = f; 

%========================================================================  
function [vertices,faces] = a_vector(ix, iy, iz, ir, hr, hl)

% ix: x coordinates of the start and end points 
% iy: y coordinates of the start and end points 
% iz: z coordinates of the start and end points 
% ir: radius of the vector line
% hr: radius of the head
% hl: head length to vector length ratio

N=36;
[x y z] = cylinder(1,N-1);
x = x(1,:).'; y = y(1,:).'; z = z(1,:).'; 

v = [[0 0 0]; x*ir y*ir z; x*ir y*ir z+(1-hl); ...
    x*hr y*hr z+(1-hl); [0 0 1]];

za = zeros(1,N-1);

f1 = [za+1; za+1; [2:N]; [3:N+1]];
f2 = [[2:N]; [3:N+1]; [3+N:2*N+1]; [2+N:2*N]];
f3 = [[2+N:2*N]; [3+N:2*N+1]; [3+2*N:3*N+1]; [2+2*N:3*N]];
f4 = [[2+2*N:3*N]; [3+2*N:3*N+1]; za+3*N+2; za+3*N+2];
faces = [f1 f2 f3 f4].';
    
% scale, rotate, and move the vector 
lenx = ix(2)-ix(1);
leny = iy(2)-iy(1);
lenz = iz(2)-iz(1);

lenxyz = (lenx^2+leny^2+lenz^2)^0.5;
lenxy = (lenx^2+leny^2)^0.5;
theta = atan2(lenxy, lenz);
phi   = atan2(leny, lenx);

v(:, 3) = v(:, 3) * lenxyz;

% rotation matrix
Rx = [1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)];
Ry = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];
Rz = [cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1];

vt = (Rz * (Ry *(v.'))).';
vt(:,1) = vt(:,1) + ix(1);
vt(:,2) = vt(:,2) + iy(1);
vt(:,3) = vt(:,3) + iz(1);
vertices = vt;

% ===============================================================
function [vertices, faces] = a_line(ix, iy, iz, ir)

% ix: x coordinates of the start and end points 
% iy: y coordinates of the start and end points 
% iz: z coordinates of the start and end points 
% ir: radius of the line

N=36;
[x y z] = cylinder(1, N-1);
x = x(1,:).'; y = y(1,:).'; z = z(1,:).'; 

v = [[0 0 0]; x*ir y*ir z; x*ir y*ir z+1; [0 0 1]];

za = ones(1,N-1);
f1 = [za; za; [2:N]; [3:N+1]];
f2 = [[2:N]; [3:N+1]; [3+N:2*N+1]; [2+N:2*N]];
f3 = [[2+N:2*N]; [3+N:2*N+1]; za+2*N+1; za+2*N+1];
faces  = [f1 f2 f3].';
    
% scale, rotate, and move the line 
lenx = ix(2)-ix(1);
leny = iy(2)-iy(1);
lenz = iz(2)-iz(1);

lenxyz = (lenx^2+leny^2+lenz^2)^0.5;
lenxy  = (lenx^2+leny^2)^0.5;
theta  = atan2(lenxy, lenz);
phi    = atan2(leny, lenx);

v(:, 3) = v(:, 3) * lenxyz;

% rotation matrix
Rx = [1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)];
Ry = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];
Rz = [cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1];

vt = (Rz * (Ry *(v.'))).';
vt(:,1) = vt(:,1) + ix(1);
vt(:,2) = vt(:,2) + iy(1);
vt(:,3) = vt(:,3) + iz(1);
vertices = vt;

% ===============================================================

function [vertices, faces] = a_sphere(cx, cy, cz, ir)
% cx: x coordinates of the sphere center
% cy: y coordinates of the sphere center
% cz: z coordinates of the sphere center
% ir: radius of the sphere

N=36;
[x,y,z] = sphere(N);
x = x * ir + cx;
y = y * ir + cy;
z = z * ir + cz;

[faces,vertices] = surf2patch(x,y,z,z);
      