function [varargout] = display_material_mesh_gui(varargin)

% A graphical user interface for plotting the material distribution

mOutputArgs     = {};       % Variable for storing output when GUI returns

handles = read_data_file;    
if isempty(handles)
    return;
end

col = [254 253 222]/255;
hMainFigure    =   figure(...       % the main GUI figure
                        'Units','characters',...
                        'MenuBar','none', ...
                        'Toolbar','none', ...
                        'HandleVisibility','callback', ...
                        'Name', 'Display Material Mesh', ...
                        'NumberTitle','off', ...
                        'Position',[5 5 100 8],...
                        'color', col);
                    
cRadioPermittivity = uicontrol('Style', 'radiobutton',...
        'Parent', hMainFigure, ...
        'String', 'relative permittivity',...
        'Value', 1,...
        'Units','characters',...
        'Position', [1 6 25 1],...
        'backgroundcolor', col, ...
        'Callback', @fRadioPermittivity);                    

cRadioPermeability = uicontrol('Style', 'radiobutton',...
        'Parent', hMainFigure, ...
        'String', 'relative permeability',...
        'Units','characters',...
        'Position', [1 4.5 25 1],...
        'backgroundcolor', col, ...
        'Callback', @fRadioPermeability);                    
    
cRadioEconductivity = uicontrol('Style', 'radiobutton',...
        'Parent', hMainFigure, ...
        'String', 'electric conductivity',...
        'Units','characters',...
        'Position', [1 3 25 1],...
        'backgroundcolor', col, ...
        'Callback', @fRadioEconductivity);                    

cRadioMconductivity = uicontrol('Style', 'radiobutton',...
        'Parent', hMainFigure, ...
        'String', 'magnetic conductivity',...
        'Units','characters',...
        'Position', [1 1.5 25 1],...
        'backgroundcolor', col, ...
        'Callback', @fRadioMconductivity);                    

cToggleXYplane = uicontrol('Style', 'toggle',...
        'Parent', hMainFigure, ...
        'String', 'xy plane',...
        'Units','characters',...
        'Position', [27 5.2 15 1.6],...
        'backgroundcolor', 'g', ...
        'value', 1,...
        'Callback', @fToggleXYplane); 
    
cToggleYZplane = uicontrol('Style', 'toggle',...
        'Parent', hMainFigure, ...
        'String', 'yz plane',...
        'Units','characters',...
        'Position', [27 3.3 15 1.6],...
        'backgroundcolor', 'g', ...
        'value', 1,...
        'Callback', @fToggleYZplane); 
    
cToggleZXplane = uicontrol('Style', 'toggle',...
        'Parent', hMainFigure, ...
        'String', 'zx plane',...
        'Units','characters',...
        'Position', [27 1.5 15 1.6],...
        'backgroundcolor', 'g', ...
        'value', 1,...
        'Callback', @fToggleZXplane); 

cEditK = uicontrol('Style', 'edit',...
        'Parent', hMainFigure, ...
        'String', 'xy plane',...
        'Units','characters',...
        'Position', [50 5.2 15 1.6],...
        'backgroundcolor', 'w', ...
        'enable','inactive', ...
        'Callback', @fEditK); 
    
cEditI = uicontrol('Style', 'edit',...
        'Parent', hMainFigure, ...
        'String', 'yz plane',...
        'Units','characters',...
        'Position', [50 3.3 15 1.6],...
        'backgroundcolor', 'w', ...
        'enable','inactive', ...
        'Callback', @fEditI); 

cEditJ = uicontrol('Style', 'edit',...
        'Parent', hMainFigure, ...
        'String', 'zx plane',...
        'Units','characters',...
        'Position', [50 1.5 15 1.6],...
        'backgroundcolor', 'w', ...
        'enable','inactive', ...
        'Callback', @fEditJ); 
    
cPushKdown = uicontrol('Style', 'pushbutton',...
        'Parent', hMainFigure, ...
        'String', '<<',...
        'Units','characters',...
        'Position', [45 5.2 5 1.6],...
        'Callback', @fPushKdown); 

cPushIdown = uicontrol('Style', 'pushbutton',...
        'Parent', hMainFigure, ...
        'String', '<<',...
        'Units','characters',...
        'Position', [45 3.3 5 1.6],...
        'Callback', @fPushIdown); 
    
cPushJdown = uicontrol('Style', 'pushbutton',...
        'Parent', hMainFigure, ...
        'String', '<<',...
        'Units','characters',...
        'Position', [45 1.5 5 1.6],...
        'Callback', @fPushJdown); 
    
cPushKup = uicontrol('Style', 'pushbutton',...
        'Parent', hMainFigure, ...
        'String', '>>',...
        'Units','characters',...
        'Position', [65 5.2 5 1.6],...
        'Callback', @fPushKup); 

cPushIup = uicontrol('Style', 'pushbutton',...
        'Parent', hMainFigure, ...
        'String', '>>',...
        'Units','characters',...
        'Position', [65 3.3 5 1.6],...
        'Callback', @fPushIup); 
    
cPushJup = uicontrol('Style', 'pushbutton',...
        'Parent', hMainFigure, ...
        'String', '>>',...
        'Units','characters',...
        'Position', [65 1.5 5 1.6],...
        'Callback', @fPushJup); 
    
cCheckEdge = uicontrol('Style', 'checkbox',...
        'Parent', hMainFigure, ...
        'String', 'edge',...
        'Units','characters',...
        'Position', [75 5.2 15 1.6],...
        'backgroundcolor', col, ...
        'Callback', @fCheckEdge); 
    
cCheckAxis = uicontrol('Style', 'checkbox',...
        'Parent', hMainFigure, ...
        'String', 'axis',...
        'Units','characters',...
        'Position', [90 5.2 15 1.6],...
        'backgroundcolor', col, ...
        'Callback', @fCheckAxis); 
    
cPushThicknessDown = uicontrol('Style', 'pushbutton',...
        'Parent', hMainFigure, ...
        'String', '<<',...
        'Units','characters',...
        'Position', [75 3.3 5 1.6],...
        'Callback', @fPushThicknessDown); 
    
cPushThicknessUp = uicontrol('Style', 'pushbutton',...
        'Parent', hMainFigure, ...
        'String', '>>',...
        'Units','characters',...
        'Position', [80 3.3 5 1.6],...
        'Callback', @fPushThicknessUp); 
    
cTextThickness = uicontrol('Style', 'text',...
        'Parent', hMainFigure, ...
        'String', 'thickness',...
        'Units','characters',...
        'backgroundcolor', col, ...
        'Position', [85 3 10 1.6]);
 
cCheckColorbar = uicontrol('Style', 'checkbox',...
        'Parent', hMainFigure, ...
        'String', 'colorbar',...
        'Units','characters',...
        'Position', [75 1.5 15 1.6],...
        'backgroundcolor', col, ...
        'Callback', @fCheckColorbar); 

set(cEditK, 'string',['k = ' num2str(round((handles.nz+1)/2)) ...
    ' / ' num2str(handles.nz+1)]);
set(cEditI, 'string',['i = ' num2str(round((handles.nx+1)/2)) ...
    ' / ' num2str(handles.nx+1)]);
set(cEditJ, 'string',['j = ' num2str(round((handles.ny+1)/2)) ...
    ' / ' num2str(handles.ny+1)]);

handles.line_thickness_factor = 6;

[handles] = initialize_plotting_arrays(handles); 
[handles] = plot_mesh(handles);
view(20,30);

mOutputArgs{1} = hMainFigure;
if nargout>0
    [varargout{1:nargout}] = mOutputArgs{:};
end

function fRadioPermittivity(hObject,eventdata)
    set(cRadioPermittivity,'value',1);
    set(cRadioPermeability,'value',0);
    set(cRadioEconductivity,'value',0);
    set(cRadioMconductivity,'value',0);
    set_maxval(handles);
    [handles] = plot_mesh(handles);    
end

function fRadioPermeability(hObject,eventdata)
    set(cRadioPermittivity,'value',0);
    set(cRadioPermeability,'value',1);
    set(cRadioEconductivity,'value',0);
    set(cRadioMconductivity,'value',0);
    set_maxval(handles);
    [handles] = plot_mesh(handles);    
end

function fRadioEconductivity(hObject,eventdata)
    set(cRadioPermittivity,'value',0);
    set(cRadioPermeability,'value',0);
    set(cRadioEconductivity,'value',1);
    set(cRadioMconductivity,'value',0);
    set_maxval(handles);
    [handles] = plot_mesh(handles);    
end

function fRadioMconductivity(hObject,eventdata)
    set(cRadioPermittivity,'value',0);
    set(cRadioPermeability,'value',0);
    set(cRadioEconductivity,'value',0);
    set(cRadioMconductivity,'value',1);
    set_maxval(handles);
    [handles] = plot_mesh(handles);
end

function fToggleXYplane(hObject,eventdata)
    val = get(hObject,'value');    
    if (val == 1) col = 'g'; else col = 'r'; end
    set(hObject,'backgroundcolor', col);
    [handles] = plot_mesh(handles);
end
function fToggleYZplane(hObject,eventdata)
    val = get(hObject,'value');    
    if (val == 1) col = 'g'; else col = 'r'; end
    set(hObject,'backgroundcolor', col);
    [handles] = plot_mesh(handles);
end
function fToggleZXplane(hObject,eventdata)
    val = get(hObject,'value');    
    if (val == 1) col = 'g'; else col = 'r'; end
    set(hObject,'backgroundcolor', col);
    [handles] = plot_mesh(handles);
end
function fEditK(hObject,eventdata)

end
function fEditI(hObject,eventdata)

end
function fEditJ(hObject,eventdata)

end
function fPushKdown(hObject,eventdata)
    
    val = get_index(get(cEditK,'string'));
    val = val-1;
    maxval = handles.nz;
    if get(cRadioPermittivity,'value')==1 || get(cRadioEconductivity,'value')==1
        maxval = maxval + 1;
    end
    if val<1
        val=1;
    end
    set(cEditK,'string',['k = ' num2str(val) ' / ' num2str(maxval)]);
    [handles] = plot_mesh(handles);

end
function fPushIdown(hObject,eventdata)
    
    val = get_index(get(cEditI,'string'));
    val = val-1;
    maxval = handles.nx;
    if get(cRadioPermittivity,'value')==1 || get(cRadioEconductivity,'value')==1
        maxval = maxval + 1;
    end
    if val<1
        val=1;
    end
    set(cEditI,'string',['i = ' num2str(val) ' / ' num2str(maxval)]);
    [handles] = plot_mesh(handles);

end
function fPushJdown(hObject,eventdata)
    val = get_index(get(cEditJ,'string'));
    val = val-1;
    maxval = handles.ny;
    if get(cRadioPermittivity,'value')==1 || get(cRadioEconductivity,'value')==1
        maxval = maxval + 1;
    end
    if val<1
        val=1;
    end
    set(cEditJ,'string',['j = ' num2str(val) ' / ' num2str(maxval)]);
    [handles] = plot_mesh(handles);
end
function fPushKup(hObject,eventdata)
    
    val = get_index(get(cEditK,'string'));
    val = val+1;
    maxval = handles.nz;
    if get(cRadioPermittivity,'value')==1 || get(cRadioEconductivity,'value')==1
        maxval = maxval + 1;
    end
    if val>maxval
        val=maxval;
    end
    set(cEditK,'string',['k = ' num2str(val) ' / ' num2str(maxval)]);
    [handles] = plot_mesh(handles);
    
end
function fPushIup(hObject,eventdata)
    
    val = get_index(get(cEditI,'string'));
    val = val+1;
    maxval = handles.nx;
    if get(cRadioPermittivity,'value')==1 || get(cRadioEconductivity,'value')==1
        maxval = maxval + 1;
    end
    if val>maxval
        val=maxval;
    end
    set(cEditI,'string',['i = ' num2str(val) ' / ' num2str(maxval)]);
    [handles] = plot_mesh(handles);

end
function fPushJup(hObject,eventdata)
    val = get_index(get(cEditJ,'string'));
    val = val+1;
    maxval = handles.ny;
    if get(cRadioPermittivity,'value')==1 || get(cRadioEconductivity,'value')==1
        maxval = maxval + 1;
    end
    if val>maxval
        val=maxval;
    end
    set(cEditJ,'string',['j = ' num2str(val) ' / ' num2str(maxval)]);
    [handles] = plot_mesh(handles);
end
function fCheckEdge(hObject,eventdata)
    [handles] = plot_mesh(handles);
end

function fCheckAxis(hObject,eventdata)
    [handles] = plot_mesh(handles);
end

function fPushThicknessDown(hObject,eventdata)
    
    handles.line_thickness_factor = handles.line_thickness_factor + 1;
    [handles] = initialize_plotting_arrays(handles);
    [handles] = plot_mesh(handles);
    guidata(hObject,handles);
    
end
function fPushThicknessUp(hObject,eventdata)
    
    handles.line_thickness_factor = handles.line_thickness_factor - 1;
    if handles.line_thickness_factor == 1
        handles.line_thickness_factor = 2;
    end
    [handles] = initialize_plotting_arrays(handles);
    [handles] = plot_mesh(handles);
    guidata(hObject,handles);
    
end

function fCheckColorbar(hObject,eventdata)
    [handles] = plot_mesh(handles);
end

function [handles] = read_data_file()

    fid = fopen('mesh_data_file.mat','r');
    if fid==-1
        msgbox('mesh_data_file does not exist');
        handles = [];
        return;
    end
    fclose(fid);
    load('mesh_data_file.mat', 'eps_r_x', 'eps_r_y', 'eps_r_z', ...
        'mu_r_x', 'mu_r_y', 'mu_r_z',...
        'sigma_e_x', 'sigma_e_y', 'sigma_e_z', ...
        'sigma_m_x', 'sigma_m_y', 'sigma_m_z', ...
        'nx','ny','nz','dx','dy','dz','fdtd_domain');
    delete('mesh_data_file.mat');

    handles.eps_r_x = eps_r_x;
    handles.eps_r_y = eps_r_y;
    handles.eps_r_z = eps_r_z;
    handles.mu_r_x = mu_r_x;
    handles.mu_r_y = mu_r_y;
    handles.mu_r_z = mu_r_z;
    handles.sigma_e_x = sigma_e_x;
    handles.sigma_e_y = sigma_e_y;
    handles.sigma_e_z = sigma_e_z;
    handles.sigma_m_x = sigma_m_x;
    handles.sigma_m_y = sigma_m_y;
    handles.sigma_m_z = sigma_m_z;
    handles.nx = nx;
    handles.ny = ny;
    handles.nz = nz;
    handles.dx = dx;
    handles.dy = dy;
    handles.dz = dz;
    handles.min_x = fdtd_domain.min_x;
    handles.min_y = fdtd_domain.min_y;
    handles.min_z = fdtd_domain.min_z;
end

function [handles] = initialize_plotting_arrays(handles)

    dx = handles.dx;
    dy = handles.dy;
    dz = handles.dz;
    nx = handles.nx;
    ny = handles.ny;
    nz = handles.nz;
    min_x = handles.min_x;
    min_y = handles.min_y;
    min_z = handles.min_z;
    
    t = handles.line_thickness_factor;
    t1 = t-1;
    
    mx = min(min(min(handles.eps_r_x)));
    my = min(min(min(handles.eps_r_y)));
    mz = min(min(min(handles.eps_r_z)));
    handles.eps_r_minval = min([mx my mz]);
    mx = max(max(max(handles.eps_r_x)));
    my = max(max(max(handles.eps_r_y)));
    mz = max(max(max(handles.eps_r_z)));
    handles.eps_r_maxval = max([mx my mz]);

    mx = min(min(min(handles.mu_r_x)));
    my = min(min(min(handles.mu_r_y)));
    mz = min(min(min(handles.mu_r_z)));
    handles.mu_r_minval = min([mx my mz]);
    mx = max(max(max(handles.mu_r_x)));
    my = max(max(max(handles.mu_r_y)));
    mz = max(max(max(handles.mu_r_z)));
    handles.mu_r_maxval = max([mx my mz]);

    mx = min(min(min(handles.sigma_e_x)));
    my = min(min(min(handles.sigma_e_y)));
    mz = min(min(min(handles.sigma_e_z)));
    handles.sigma_e_minval = min([mx my mz]);
    mx = max(max(max(handles.sigma_e_x)));
    my = max(max(max(handles.sigma_e_y)));
    mz = max(max(max(handles.sigma_e_z)));
    handles.sigma_e_maxval = max([mx my mz]);

    mx = min(min(min(handles.sigma_m_x)));
    my = min(min(min(handles.sigma_m_y)));
    mz = min(min(min(handles.sigma_m_z)));
    handles.sigma_m_minval = min([mx my mz]);
    mx = max(max(max(handles.sigma_m_x)));
    my = max(max(max(handles.sigma_m_y)));
    mz = max(max(max(handles.sigma_m_z)));
    handles.sigma_m_maxval = max([mx my mz]);

    %eps_r xy plane x component
    npatches = nx*(ny+1);
    fac = 1:npatches*6;
    fac = reshape(fac,6,[]).';
    handles.faces_xy_x = fac;
    
    x = meshgrid(0:nx-1,1:ny+1)*dx;
    x = reshape(x,1,[]).';
    x = [x x+dx/t x+t1*dx/t x+dx x+t1*dx/t x+dx/t];
    x = reshape(x.',1,[]).';
    handles.vertices_xy_x_x = x + min_x;
    
    y = (meshgrid(0:ny,1:nx)*dy).';
    y = reshape(y,1,[]).';
    y = [y y-dy/t y-dy/t y y+dy/t y+dy/t];
    y = reshape(y.',1,[]).';
    handles.vertices_xy_x_y = y + min_y;

    %eps_r xy plane y component
    npatches = ny*(nx+1);
    fac = 1:npatches*6;
    fac = reshape(fac,6,[]).';
    handles.faces_xy_y = fac;
    
    y = meshgrid(0:ny-1,1:nx+1)*dy;
    y = reshape(y,1,[]).';
    y = [y y+dy/t y+t1*dy/t y+dy y+t1*dy/t y+dy/t];
    y = reshape(y.',1,[]).';
    handles.vertices_xy_y_y = y + min_y; 
    
    x = (meshgrid(0:nx,1:ny)*dx).';
    x = reshape(x,1,[]).';
    x = [x x-dx/t x-dx/t x x+dx/t x+dx/t];
    x = reshape(x.',1,[]).';
    handles.vertices_xy_y_x = x + min_x;

    %eps_r yz plane y component
    npatches = ny*(nz+1);
    fac = 1:npatches*6;
    fac = reshape(fac,6,[]).';
    handles.faces_yz_y = fac;
    
    y = meshgrid(0:ny-1,1:nz+1)*dy;
    y = reshape(y,1,[]).';
    y = [y y+dy/t y+t1*dy/t y+dy y+t1*dy/t y+dy/t];
    y = reshape(y.',1,[]).';
    handles.vertices_yz_y_y = y + min_y;
    
    z = (meshgrid(0:nz,1:ny)*dz).';
    z = reshape(z,1,[]).';
    z = [z z-dz/t z-dz/t z z+dz/t z+dz/t];
    z = reshape(z.',1,[]).';
    handles.vertices_yz_y_z = z + min_z;

    %eps_r yz plane z component
    npatches = nz*(ny+1);
    fac = 1:npatches*6;
    fac = reshape(fac,6,[]).';
    handles.faces_yz_z = fac;
    
    z = meshgrid(0:nz-1,1:ny+1)*dz;
    z = reshape(z,1,[]).';
    z = [z z+dz/t z+t1*dz/t z+dz z+t1*dz/t z+dz/t];
    z = reshape(z.',1,[]).';
    handles.vertices_yz_z_z = z + min_z; 
    
    y = (meshgrid(0:ny,1:nz)*dy).';
    y = reshape(y,1,[]).';
    y = [y y-dy/t y-dy/t y y+dy/t y+dy/t];
    y = reshape(y.',1,[]).';
    handles.vertices_yz_z_y = y + min_y;

    %eps_r zx plane z component
    npatches = nz*(nx+1);
    fac = 1:npatches*6;
    fac = reshape(fac,6,[]).';
    handles.faces_zx_z = fac;
    
    z = meshgrid(0:nz-1,1:nx+1)*dz;
    z = reshape(z,1,[]).';
    z = [z z+dz/t z+t1*dz/t z+dz z+t1*dz/t z+dz/t];
    z = reshape(z.',1,[]).';
    handles.vertices_zx_z_z = z + min_z;
    
    x = (meshgrid(0:nx,1:nz)*dx).';
    x = reshape(x,1,[]).';
    x = [x x-dx/t x-dx/t x x+dx/t x+dx/t];
    x = reshape(x.',1,[]).';
    handles.vertices_zx_z_x = x + min_x;

    %eps_r zx plane x component
    npatches = nx*(nz+1);
    fac = 1:npatches*6;
    fac = reshape(fac,6,[]).';
    handles.faces_zx_x = fac;
    
    x = meshgrid(0:nx-1,1:nz+1)*dx;
    x = reshape(x,1,[]).';
    x = [x x+dx/t x+t1*dx/t x+dx x+t1*dx/t x+dx/t];
    x = reshape(x.',1,[]).';
    handles.vertices_zx_x_x = x + min_x; 
    
    z = (meshgrid(0:nz,1:nx)*dz).';
    z = reshape(z,1,[]).';
    z = [z z-dz/t z-dz/t z z+dz/t z+dz/t];
    z = reshape(z.',1,[]).';
    handles.vertices_zx_x_z = z + min_z;

    %mu_r xy plane x component
    npatches = (nx+1)*ny;
    fac = 1:npatches*6;
    fac = reshape(fac,6,[]).';
    handles.mu_faces_xy_x = fac;
    
    x = meshgrid(0:nx,1:ny)*dx;
    x = reshape(x,1,[]).';
    x = [x x+dx/t x+t1*dx/t x+dx x+t1*dx/t x+dx/t];
    x = reshape(x.',1,[]).';
    handles.mu_vertices_xy_x_x = x-dx/2 + min_x;
    
    y = (meshgrid(0:ny-1,1:nx+1)*dy).';
    y = reshape(y,1,[]).';
    y = [y y-dy/t y-dy/t y y+dy/t y+dy/t];
    y = reshape(y.',1,[]).';
    handles.mu_vertices_xy_x_y = y+dy/2 + min_y;

    %mu_r xy plane y component
    npatches = (ny+1)*nx;
    fac = 1:npatches*6;
    fac = reshape(fac,6,[]).';
    handles.mu_faces_xy_y = fac;
    
    y = meshgrid(0:ny,1:nx)*dy;
    y = reshape(y,1,[]).';
    y = [y y+dy/t y+t1*dy/t y+dy y+t1*dy/t y+dy/t];
    y = reshape(y.',1,[]).';
    handles.mu_vertices_xy_y_y = y-dy/2 + min_y;
    
    x = (meshgrid(0:nx-1,1:ny+1)*dx).';
    x = reshape(x,1,[]).';
    x = [x x-dx/t x-dx/t x x+dx/t x+dx/t];
    x = reshape(x.',1,[]).';
    handles.mu_vertices_xy_y_x = x+dx/2 + min_x;

    %mu_r yz plane y component
    npatches = (ny+1)*nz;
    fac = 1:npatches*6;
    fac = reshape(fac,6,[]).';
    handles.mu_faces_yz_y = fac;
    
    y = meshgrid(0:ny,1:nz)*dy;
    y = reshape(y,1,[]).';
    y = [y y+dy/t y+t1*dy/t y+dy y+t1*dy/t y+dy/t];
    y = reshape(y.',1,[]).';
    handles.mu_vertices_yz_y_y = y-dy/2 + min_y;
    
    z = (meshgrid(0:nz-1,1:ny+1)*dz).';
    z = reshape(z,1,[]).';
    z = [z z-dz/t z-dz/t z z+dz/t z+dz/t];
    z = reshape(z.',1,[]).';
    handles.mu_vertices_yz_y_z = z+dz/2 + min_z;

    %mu_r yz plane z component
    npatches = (nz+1)*ny;
    fac = 1:npatches*6;
    fac = reshape(fac,6,[]).';
    handles.mu_faces_yz_z = fac;
    
    z = meshgrid(0:nz,1:ny)*dz;
    z = reshape(z,1,[]).';
    z = [z z+dz/t z+t1*dz/t z+dz z+t1*dz/t z+dz/t];
    z = reshape(z.',1,[]).';
    handles.mu_vertices_yz_z_z = z-dz/2 + min_z;
    
    y = (meshgrid(0:ny-1,1:nz+1)*dy).';
    y = reshape(y,1,[]).';
    y = [y y-dy/t y-dy/t y y+dy/t y+dy/t];
    y = reshape(y.',1,[]).';
    handles.mu_vertices_yz_z_y = y+dy/2 + min_y;

    %mu_r zx plane z component
    npatches = (nz+1)*nx;
    fac = 1:npatches*6;
    fac = reshape(fac,6,[]).';
    handles.mu_faces_zx_z = fac;
    
    z = meshgrid(0:nz,1:nx)*dz;
    z = reshape(z,1,[]).';
    z = [z z+dz/t z+t1*dz/t z+dz z+t1*dz/t z+dz/t];
    z = reshape(z.',1,[]).';
    handles.mu_vertices_zx_z_z = z-dz/2 + min_z;
    
    x = (meshgrid(0:nx-1,1:nz+1)*dx).';
    x = reshape(x,1,[]).';
    x = [x x-dx/t x-dx/t x x+dx/t x+dx/t];
    x = reshape(x.',1,[]).';
    handles.mu_vertices_zx_z_x = x+dx/2 + min_x;

    %mu_r zx plane x component
    npatches = (nx+1)*nz;
    fac = 1:npatches*6;
    fac = reshape(fac,6,[]).';
    handles.mu_faces_zx_x = fac;
    
    x = meshgrid(0:nx,1:nz)*dx;
    x = reshape(x,1,[]).';
    x = [x x+dx/t x+t1*dx/t x+dx x+t1*dx/t x+dx/t];
    x = reshape(x.',1,[]).';
    handles.mu_vertices_zx_x_x = x-dx/2 + min_x;
    
    z = (meshgrid(0:nz-1,1:nx+1)*dz).';
    z = reshape(z,1,[]).';
    z = [z z-dz/t z-dz/t z z+dz/t z+dz/t];
    z = reshape(z.',1,[]).';
    handles.mu_vertices_zx_x_z = z+dz/2 + min_z;
end

function [handles] = plot_mesh(handles)
    f = figure(111111);
    set(f,'name','Material Mesh');
    set(f,'numbertitle','off');
    cameratoolbar(f);
    po = findobj(gca,'type','patch');
    delete(po);
    daspect([1 1 1]);

    if get(cRadioPermittivity,'value')==1 
        if get(cToggleXYplane,'value')==1 
            [handles] = plot_eps_r_xy(handles);
        end
        if get(cToggleYZplane,'value')==1 
            [handles] = plot_eps_r_yz(handles);
        end
        if get(cToggleZXplane,'value')==1 
            [handles] = plot_eps_r_zx(handles);
        end
    end
    if get(cRadioPermeability,'value')==1 
        if get(cToggleXYplane,'value')==1 
            [handles] = plot_mu_r_xy(handles);
        end
        if get(cToggleYZplane,'value')==1 
            [handles] = plot_mu_r_yz(handles);
        end
        if get(cToggleZXplane,'value')==1 
            [handles] = plot_mu_r_zx(handles);
        end
    end
    if get(cRadioEconductivity,'value')==1 
        if get(cToggleXYplane,'value')==1 
            [handles] = plot_sigma_e_xy(handles);
        end
        if get(cToggleYZplane,'value')==1 
            [handles] = plot_sigma_e_yz(handles);
        end
        if get(cToggleZXplane,'value')==1 
            [handles] = plot_sigma_e_zx(handles);
        end
    end
    if get(cRadioMconductivity,'value')==1 
        if get(cToggleXYplane,'value')==1 
            [handles] = plot_sigma_m_xy(handles);
        end
        if get(cToggleYZplane,'value')==1 
            [handles] = plot_sigma_m_yz(handles);
        end
        if get(cToggleZXplane,'value')==1 
            [handles] = plot_sigma_m_zx(handles);
        end
    end

    if (get(cCheckEdge,'value') == 0)
        po = findobj(gca,'type','patch');
        set(po,'EdgeColor','none');
    end
    
    if (get(cCheckColorbar,'value') == 1)
        colorbar;
    else
        colorbar off;
    end    
    
    if (get(cCheckAxis,'value') == 1)
        axis on;
    else
        axis off;
    end    
end

function [handles] = plot_eps_r_xy(handles)

    val = get_index(get(cEditK,'string'));
    colormap(jet);
    caxis([handles.eps_r_minval*0.9 handles.eps_r_maxval*1.1]);

    % plot x components;
    pdata = handles.eps_r_x(:,:,val);
    zz = handles.min_z+(val-1)*handles.dz;
    x = handles.vertices_xy_x_x;
    y = handles.vertices_xy_x_y;
    z = y*0+zz;
    pdata = reshape(pdata.',1,[]);
    patch('faces',handles.faces_xy_x,'vertices',[x y z],'facecolor','flat','cdata',pdata);

    % plot y components;
    pdata = handles.eps_r_y(:,:,val);
    zz = handles.min_z+(val-1)*handles.dz;
    x = handles.vertices_xy_y_x;
    y = handles.vertices_xy_y_y;
    z = x*0+zz;
    pdata = reshape(pdata,1,[]);
    patch('faces',handles.faces_xy_y,'vertices',[x y z],'facecolor','flat','cdata',pdata);

end

function [handles] = plot_eps_r_yz(handles)
    
    val = get_index(get(cEditI,'string'));
    colormap(jet);
    caxis([handles.eps_r_minval*0.9 handles.eps_r_maxval*1.1]);

    % plot y components;
    pdata = squeeze(handles.eps_r_y(val,:,:));
    xx = handles.min_x+(val-1)*handles.dx;
    y = handles.vertices_yz_y_y;
    z = handles.vertices_yz_y_z;
    x = z*0+xx;
    pdata = reshape(pdata.',1,[]);
    patch('faces',handles.faces_yz_y,'vertices',[x y z],'facecolor','flat','cdata',pdata);

    % plot z components;
    pdata = squeeze(handles.eps_r_z(val,:,:));
    xx = handles.min_x+(val-1)*handles.dx;
    y = handles.vertices_yz_z_y;
    z = handles.vertices_yz_z_z;
    x = y*0+xx;
    pdata = reshape(pdata,1,[]);
    patch('faces',handles.faces_yz_z,'vertices',[x y z],'facecolor','flat','cdata',pdata);

end

function [handles] = plot_eps_r_zx(handles)

    val = get_index(get(cEditJ,'string'));
    colormap(jet);
    caxis([handles.eps_r_minval*0.9 handles.eps_r_maxval*1.1]);

    % plot z components;
    pdata = squeeze(handles.eps_r_z(:,val,:));
    yy = handles.min_y+(val-1)*handles.dy;
    z = handles.vertices_zx_z_z;
    x = handles.vertices_zx_z_x;
    y = x*0+yy;
    pdata = reshape(pdata,1,[]);
    patch('faces',handles.faces_zx_z,'vertices',[x y z],'facecolor','flat','cdata',pdata);

    % plot x components;
    pdata = squeeze(handles.eps_r_x(:,val,:));
    yy = handles.min_y+(val-1)*handles.dy;
    z = handles.vertices_zx_x_z;
    x = handles.vertices_zx_x_x;
    y = z*0+yy;
    pdata = reshape(pdata.',1,[]);
    patch('faces',handles.faces_zx_x,'vertices',[x y z],'facecolor','flat','cdata',pdata);
end

function [handles] = plot_sigma_e_xy(handles)

    val = get_index(get(cEditK,'string'));
    colormap(jet);
    caxis([handles.sigma_e_minval*0.9 handles.sigma_e_maxval*1.1]);

    % plot x components;
    pdata = handles.sigma_e_x(:,:,val);
    zz = handles.min_z+(val-1)*handles.dz;
    x = handles.vertices_xy_x_x;
    y = handles.vertices_xy_x_y;
    z = y*0+zz;
    pdata = reshape(pdata.',1,[]);
    patch('faces',handles.faces_xy_x,'vertices',[x y z],'facecolor','flat','cdata',pdata);

    % plot y components;
    pdata = handles.sigma_e_y(:,:,val);
    zz = handles.min_z+(val-1)*handles.dz;
    x = handles.vertices_xy_y_x;
    y = handles.vertices_xy_y_y;
    z = x*0+zz;
    pdata = reshape(pdata,1,[]);
    patch('faces',handles.faces_xy_y,'vertices',[x y z],'facecolor','flat','cdata',pdata);
end

function [handles] = plot_sigma_e_yz(handles)
    
    val = get_index(get(cEditI,'string'));
    colormap(jet);
    caxis([handles.sigma_e_minval*0.9 handles.sigma_e_maxval*1.1]);

    % plot y components;
    pdata = squeeze(handles.sigma_e_y(val,:,:));
    xx = handles.min_x+(val-1)*handles.dx;
    y = handles.vertices_yz_y_y;
    z = handles.vertices_yz_y_z;
    x = z*0+xx;
    pdata = reshape(pdata.',1,[]);
    patch('faces',handles.faces_yz_y,'vertices',[x y z],'facecolor','flat','cdata',pdata);

    % plot z components;
    pdata = squeeze(handles.sigma_e_z(val,:,:));
    xx = handles.min_x+(val-1)*handles.dx;
    y = handles.vertices_yz_z_y;
    z = handles.vertices_yz_z_z;
    x = y*0+xx;
    pdata = reshape(pdata,1,[]);
    patch('faces',handles.faces_yz_z,'vertices',[x y z],'facecolor','flat','cdata',pdata);

end

function [handles] = plot_sigma_e_zx(handles)

    val = get_index(get(cEditJ,'string'));
    colormap(jet);
    caxis([handles.sigma_e_minval*0.9 handles.sigma_e_maxval*1.1]);

    % plot z components;
    pdata = squeeze(handles.sigma_e_z(:,val,:));
    yy = handles.min_y+(val-1)*handles.dy;
    z = handles.vertices_zx_z_z;
    x = handles.vertices_zx_z_x;
    y = x*0+yy;
    pdata = reshape(pdata,1,[]);
    patch('faces',handles.faces_zx_z,'vertices',[x y z],'facecolor','flat','cdata',pdata);

    % plot x components;
    pdata = squeeze(handles.sigma_e_x(:,val,:));
    yy = handles.min_y+(val-1)*handles.dy;
    z = handles.vertices_zx_x_z;
    x = handles.vertices_zx_x_x;
    y = z*0+yy;
    pdata = reshape(pdata.',1,[]);
    patch('faces',handles.faces_zx_x,'vertices',[x y z],'facecolor','flat','cdata',pdata);

end

function [handles] = plot_mu_r_xy(handles)

    val = get_index(get(cEditK,'string'));
    colormap(jet);
    caxis([handles.mu_r_minval*0.9 handles.mu_r_maxval*1.1]);

    % plot x components;
    pdata = handles.mu_r_x(:,:,val);
    zz = handles.min_z+(val-0.5)*handles.dz;
    x = handles.mu_vertices_xy_x_x;
    y = handles.mu_vertices_xy_x_y;
    z = y*0+zz;
    pdata = reshape(pdata.',1,[]);
    patch('faces',handles.mu_faces_xy_x,'vertices',[x y z],'facecolor','flat','cdata',pdata);

    % plot y components;
    pdata = handles.mu_r_y(:,:,val);
    zz = handles.min_z+(val-0.5)*handles.dz;
    x = handles.mu_vertices_xy_y_x;
    y = handles.mu_vertices_xy_y_y;
    z = x*0+zz;
    pdata = reshape(pdata,1,[]);
    patch('faces',handles.mu_faces_xy_y,'vertices',[x y z],'facecolor','flat','cdata',pdata);
end

function [handles] = plot_mu_r_yz(handles)
    val = get_index(get(cEditI,'string'));
    colormap(jet);
    caxis([handles.mu_r_minval*0.9 handles.mu_r_maxval*1.1]);

    % plot y components;
    pdata = squeeze(handles.mu_r_y(val,:,:));
    xx = handles.min_x+(val-0.5)*handles.dx;
    y = handles.mu_vertices_yz_y_y;
    z = handles.mu_vertices_yz_y_z;
    x = z*0+xx;
    pdata = reshape(pdata.',1,[]);
    patch('faces',handles.mu_faces_yz_y,'vertices',[x y z],'facecolor','flat','cdata',pdata);

    % plot z components;
    pdata = squeeze(handles.mu_r_z(val,:,:));
    xx = handles.min_x+(val-0.5)*handles.dx;
    y = handles.mu_vertices_yz_z_y;
    z = handles.mu_vertices_yz_z_z;
    x = y*0+xx;
    pdata = reshape(pdata,1,[]);
    patch('faces',handles.mu_faces_yz_z,'vertices',[x y z],'facecolor','flat','cdata',pdata);
end

function [handles] = plot_mu_r_zx(handles)

    val = get_index(get(cEditJ,'string'));
    colormap(jet);
    caxis([handles.mu_r_minval*0.9 handles.mu_r_maxval*1.1]);

    % plot z components;
    pdata = squeeze(handles.mu_r_z(:,val,:));
    yy = handles.min_y+(val-0.5)*handles.dy;
    z = handles.mu_vertices_zx_z_z;
    x = handles.mu_vertices_zx_z_x;
    y = x*0+yy;
    pdata = reshape(pdata,1,[]);
    patch('faces',handles.mu_faces_zx_z,'vertices',[x y z],'facecolor','flat','cdata',pdata);

    % plot x components;
    pdata = squeeze(handles.mu_r_x(:,val,:));
    yy = handles.min_y+(val-0.5)*handles.dy;
    z = handles.mu_vertices_zx_x_z;
    x = handles.mu_vertices_zx_x_x;
    y = z*0+yy;
    pdata = reshape(pdata.',1,[]);
    patch('faces',handles.mu_faces_zx_x,'vertices',[x y z],'facecolor','flat','cdata',pdata);
end

function [handles] = plot_sigma_m_xy(handles)

    val = get_index(get(cEditK,'string'));
    colormap(jet);
    caxis([handles.sigma_m_minval*0.9 handles.sigma_m_maxval*1.1]);

    % plot x components;
    pdata = handles.sigma_m_x(:,:,val);
    zz = handles.min_z+(val-0.5)*handles.dz;
    x = handles.mu_vertices_xy_x_x;
    y = handles.mu_vertices_xy_x_y;
    z = y*0+zz;
    pdata = reshape(pdata.',1,[]);
    patch('faces',handles.mu_faces_xy_x,'vertices',[x y z],'facecolor','flat','cdata',pdata);

    % plot y components;
    pdata = handles.sigma_m_y(:,:,val);
    zz = handles.min_z+(val-0.5)*handles.dz;
    x = handles.mu_vertices_xy_y_x;
    y = handles.mu_vertices_xy_y_y;
    z = x*0+zz;
    pdata = reshape(pdata,1,[]);
    patch('faces',handles.mu_faces_xy_y,'vertices',[x y z],'facecolor','flat','cdata',pdata);
end

function [handles] = plot_sigma_m_yz(handles)
    
    val = get_index(get(cEditI,'string'));
    colormap(jet);
    caxis([handles.sigma_m_minval*0.9 handles.sigma_m_maxval*1.1]);

    % plot y components;
    pdata = squeeze(handles.sigma_m_y(val,:,:));
    xx = handles.min_x+(val-0.5)*handles.dx;
    y = handles.mu_vertices_yz_y_y;
    z = handles.mu_vertices_yz_y_z;
    x = z*0+xx;
    pdata = reshape(pdata.',1,[]);
    patch('faces',handles.mu_faces_yz_y,'vertices',[x y z],'facecolor','flat','cdata',pdata);

    % plot z components;
    pdata = squeeze(handles.sigma_m_z(val,:,:));
    xx = handles.min_x+(val-0.5)*handles.dx;
    y = handles.mu_vertices_yz_z_y;
    z = handles.mu_vertices_yz_z_z;
    x = y*0+xx;
    pdata = reshape(pdata,1,[]);
    patch('faces',handles.mu_faces_yz_z,'vertices',[x y z],'facecolor','flat','cdata',pdata);
end

function [handles] = plot_sigma_m_zx(handles)

    val = get_index(get(cEditJ,'string'));
    colormap(jet);
    caxis([handles.sigma_m_minval*0.9 handles.sigma_m_maxval*1.1]);

    % plot z components;
    pdata = squeeze(handles.sigma_m_z(:,val,:));
    yy = handles.min_y+(val-0.5)*handles.dy;
    z = handles.mu_vertices_zx_z_z;
    x = handles.mu_vertices_zx_z_x;
    y = x*0+yy;
    pdata = reshape(pdata,1,[]);
    patch('faces',handles.mu_faces_zx_z,'vertices',[x y z],'facecolor','flat','cdata',pdata);

    % plot x components;
    pdata = squeeze(handles.sigma_m_x(:,val,:));
    yy = handles.min_y+(val-0.5)*handles.dy;
    z = handles.mu_vertices_zx_x_z;
    x = handles.mu_vertices_zx_x_x;
    y = z*0+yy;
    pdata = reshape(pdata.',1,[]);
    patch('faces',handles.mu_faces_zx_x,'vertices',[x y z],'facecolor','flat','cdata',pdata);
end

function [val] = get_index(str)
    [token, remain] = strtok(str,'=');
    remain = strrep(remain, '=', ' ');
    [token, remain] = strtok(remain,'/');
    val = str2num(token);
end
function [] = set_maxval(handles)
    
    val = get_index(get(cEditJ,'string'));
    maxval = handles.ny;
    if get(cRadioPermittivity,'value')==1 || get(cRadioEconductivity,'value')==1
        maxval = maxval + 1;
    end
    if val>maxval
        val=maxval;
    end
    set(cEditJ,'string',['j = ' num2str(val) ' / ' num2str(maxval)]);

    val = get_index(get(cEditI,'string'));
    maxval = handles.nx;
    if get(cRadioPermittivity,'value')==1 || get(cRadioEconductivity,'value')==1
        maxval = maxval + 1;
    end
    if val>maxval
        val=maxval;
    end
    set(cEditI,'string',['i = ' num2str(val) ' / ' num2str(maxval)]);

    val = get_index(get(cEditK,'string'));
    maxval = handles.nz;
    if get(cRadioPermittivity,'value')==1 || get(cRadioEconductivity,'value')==1
        maxval = maxval + 1;
    end
    if val>maxval
        val=maxval;
    end
    set(cEditK,'string',['k = ' num2str(val) ' / ' num2str(maxval)]);
end

end
