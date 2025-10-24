
% display voltage sources
if exist('voltage_sources','var')
    for i=1:size(voltage_sources,2)
        obj = voltage_sources(i);
        obj.obj_type = 'voltage_source';
        obj.obj_id = ['Vs_' num2str(i)];
        object_drawing_functions(obj);
    end
end

% display current sources
if exist('current_sources','var')
    for i=1:size(current_sources,2)
        obj = current_sources(i);
        obj.obj_type = 'current_source';
        obj.obj_id = ['Is_' num2str(i)];
        object_drawing_functions(obj);
    end
end

% display resistors
if exist('resistors','var')
    for i=1:size(resistors,2)
        obj = resistors(i);
        obj.obj_type = 'resistor';
        obj.obj_id = ['R_' num2str(i)];
        object_drawing_functions(obj);
    end
end

% display capacitors
if exist('capacitors','var')
    for i=1:size(capacitors,2)
        obj = capacitors(i);
        obj.obj_type = 'capacitor';
        obj.obj_id = ['C_' num2str(i)];
        object_drawing_functions(obj);
    end
end

% display inductors
if exist('inductors','var')
    for i=1:size(inductors,2)
        obj = inductors(i);
        obj.obj_type = 'inductor';
        obj.obj_id = ['L_' num2str(i)];
        object_drawing_functions(obj);
    end
end

% display diodes
if exist('diodes','var')
    for i=1:size(diodes,2)
        obj = diodes(i);
        obj.obj_type = 'diode';
        obj.obj_id = ['D_' num2str(i)];
        object_drawing_functions(obj);
    end
end

% display sampled voltages
if exist('sampled_voltages','var')
    for i=1:size(sampled_voltages,2)
        obj = sampled_voltages(i);
        obj.obj_type = 'sampled_voltage';
        obj.obj_id = ['SV_' num2str(i)];
        object_drawing_functions(obj);
    end
end

% display sampled currents
if exist('sampled_currents','var')
    for i=1:size(sampled_currents,2)
        obj = sampled_currents(i);
        obj.obj_type = 'sampled_current';
        obj.obj_id = ['SC_' num2str(i)];
        object_drawing_functions(obj);
    end
end

% display thin wires
if exist('thin_wires','var')
    for i=1:size(thin_wires,2)
        obj = thin_wires(i);
        obj.obj_type = 'thin_wire';
        obj.obj_id = ['TW_' num2str(i)];
        object_drawing_functions(obj);
    end
end

% display sampled electric fields
if exist('sampled_electric_fields','var')
    for i=1:size(sampled_electric_fields,2)
        obj = sampled_electric_fields(i);
        obj.obj_type = 'sampled_electric_field';
        obj.obj_id = ['sef_' num2str(i)];
        obj.dx = dx; obj.dy = dy; obj.dz = dz;
        object_drawing_functions(obj);
    end
end

% display sampled magnetic fields
if exist('sampled_magnetic_fields','var')
    for i=1:size(sampled_magnetic_fields,2)
        obj = sampled_magnetic_fields(i);
        obj.obj_type = 'sampled_magnetic_field';
        obj.obj_id = ['smf_' num2str(i)];
        obj.dx = dx; obj.dy = dy; obj.dz = dz;
        object_drawing_functions(obj);
    end
end

% draw axes outside the domain
if problem_space_display.axis_outside_domain
    obj.obj_type = 'axis_outside_domain';
    obj.lx = fdtd_domain.min_x;
    obj.ly = fdtd_domain.min_y;
    obj.lz = fdtd_domain.min_z;
    obj.dx = dx; obj.dy = dy; obj.dz = dz;
    object_drawing_functions(obj);
end

% draw axes at the origin
if problem_space_display.axis_at_origin
    obj.obj_type = 'axis_at_origin';
    obj.dx = dx; obj.dy = dy; obj.dz = dz;
    object_drawing_functions(obj);
end

% draw boundaries
if problem_space_display.outer_boundaries
    lx = fdtd_domain.min_x;
    ly = fdtd_domain.min_y;
    lz = fdtd_domain.min_z;
    ux = nx*dx+lx;
    uy = ny*dy+ly;
    uz = nz*dz+lz;
    line([lx lx], [ly ly], [lz uz], 'linewidth', 2, 'color', 'b');
    line([lx lx], [uy uy], [lz uz], 'linewidth', 2, 'color', 'b');
    line([ux ux], [ly ly], [lz uz], 'linewidth', 2, 'color', 'b');
    line([ux ux], [uy uy], [lz uz], 'linewidth', 2, 'color', 'b');
    line([lx lx], [ly uy], [lz lz], 'linewidth', 2, 'color', 'b');
    line([lx lx], [ly uy], [uz uz], 'linewidth', 2, 'color', 'b');
    line([ux ux], [ly uy], [lz lz], 'linewidth', 2, 'color', 'b');
    line([ux ux], [ly uy], [uz uz], 'linewidth', 2, 'color', 'b');
    line([lx ux], [ly ly], [lz lz], 'linewidth', 2, 'color', 'b');
    line([lx ux], [ly ly], [uz uz], 'linewidth', 2, 'color', 'b');
    line([lx ux], [uy uy], [lz lz], 'linewidth', 2, 'color', 'b');
    line([lx ux], [uy uy], [uz uz], 'linewidth', 2, 'color', 'b');
end

% draw cpml boundaries
if problem_space_display.cpml_boundaries
    lx = fdtd_domain.min_x;
    ly = fdtd_domain.min_y;
    lz = fdtd_domain.min_z;
    ux = nx*dx+lx;
    uy = ny*dy+ly;
    uz = nz*dz+lz;
    if strcmp(boundary.type_xn, 'cpml')
        lx = lx + abs(dx * boundary.cpml_number_of_cells_xn);
    end
    if strcmp(boundary.type_yn, 'cpml')
        ly = ly + abs(dy * boundary.cpml_number_of_cells_yn);
    end
    if strcmp(boundary.type_zn, 'cpml')
        lz = lz + abs(dz * boundary.cpml_number_of_cells_zn);
    end
    if strcmp(boundary.type_xp, 'cpml')
        ux = ux - abs(dx * boundary.cpml_number_of_cells_xp);
    end
    if strcmp(boundary.type_yp, 'cpml')
        uy = uy - abs(dy * boundary.cpml_number_of_cells_yp);
    end
    if strcmp(boundary.type_zp, 'cpml')
        uz = uz - abs(dz * boundary.cpml_number_of_cells_zp);
    end
    line([lx lx], [ly ly], [lz uz], 'linewidth', 2, ...
        'color', 'r', 'linestyle', '--');
    line([lx lx], [uy uy], [lz uz], 'linewidth', 2, ...
        'color', 'r', 'linestyle', '--');
    line([ux ux], [ly ly], [lz uz], 'linewidth', 2, 'color', ...
        'r', 'linestyle', '--');
    line([ux ux], [uy uy], [lz uz], 'linewidth', 2, 'color', ...
        'r', 'linestyle', '--');
    line([lx lx], [ly uy], [lz lz], 'linewidth', 2, 'color', ...
        'r', 'linestyle', '--');
    line([lx lx], [ly uy], [uz uz], 'linewidth', 2, 'color', 'r', ...
        'linestyle', '--');
    line([ux ux], [ly uy], [lz lz], 'linewidth', 2, 'color', 'r', ...
        'linestyle', '--');
    line([ux ux], [ly uy], [uz uz], 'linewidth', 2, 'color', 'r', ...
        'linestyle', '--');
    line([lx ux], [ly ly], [lz lz], 'linewidth', 2, 'color', 'r', ...
        'linestyle', '--');
    line([lx ux], [ly ly], [uz uz], 'linewidth', 2, 'color', 'r', ...
        'linestyle', '--');
    line([lx ux], [uy uy], [lz lz], 'linewidth', 2, 'color', 'r', ...
        'linestyle', '--');
    line([lx ux], [uy uy], [uz uz], 'linewidth', 2, 'color', 'r', ...
        'linestyle', '--');
end

% draw grid
lx = fdtd_domain.min_x;
ly = fdtd_domain.min_y;
lz = fdtd_domain.min_z;
ux = nx*dx+lx;
uy = ny*dy+ly;
uz = nz*dz+lz;

cl = [0.5 0.5 0.5];
if problem_space_display.grid_zn
    for i=0:ny
         line([lx ux],[ly+i*dy ly+i*dy], [lz lz],'color',cl,'linewidth',0.5);
    end
    for i=0:nx
         line([lx+i*dx lx+i*dx], [ly uy], [lz lz],'color',cl,'linewidth',0.5);
    end
end
if problem_space_display.grid_zp
    for i=0:ny
         line([lx ux],[ly+i*dy ly+i*dy], [uz uz],'color',cl,'linewidth',0.5);
    end
    for i=0:nx
         line([lx+i*dx lx+i*dx], [ly uy], [uz uz],'color',cl,'linewidth',0.5);
    end
end
if problem_space_display.grid_xn
    for i=0:nz
         line([lx lx],[ly uy],[lz+i*dz lz+i*dz],'color',cl,'linewidth',0.5);
    end
    for i=0:ny
         line([lx lx],[ly+i*dy ly+i*dy], [lz uz],'color',cl,'linewidth',0.5);
    end
end    
if problem_space_display.grid_xp
    for i=0:nz
         line([ux ux],[ly uy],[lz+i*dz lz+i*dz],'color',cl,'linewidth',0.5);
    end
    for i=0:ny
         line([ux ux],[ly+i*dy ly+i*dy], [lz uz],'color',cl,'linewidth',0.5);
    end
end    
if problem_space_display.grid_yp
    for i=0:nx
         line([lx+i*dx lx+i*dx],[uy uy],[lz uz],'color',cl,'linewidth',0.5);
    end
    for i=0:nz
         line([lx ux],[uy uy],[lz+i*dz lz+i*dz],'color',cl,'linewidth',0.5);
    end
end
if problem_space_display.grid_yn
    for i=0:nx
         line([lx+i*dx lx+i*dx],[ly ly],[lz uz],'color',cl,'linewidth',0.5);
    end
    for i=0:nz
         line([lx ux],[ly ly],[lz+i*dz lz+i*dz],'color',cl,'linewidth',0.5);
    end
end

% delete the labels
if ~problem_space_display.labels
    lbls = findobj(gcf,'type','text','-and','-not','edgecolor','none');
    delete(lbls);
end
