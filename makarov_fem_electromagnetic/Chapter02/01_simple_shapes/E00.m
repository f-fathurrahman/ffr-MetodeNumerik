function fig_hdl = E00
%   SYNTAX
%   E00
%   DESCRIPTION
%   This module outputs uniform and non-uniform triangular meshes for a plate,
%   circle, brick, sphere, or a cylinder arbitrarily oriented in space into
%   the *.mat file mesh.mat. The mesh includes an array of nodes P(:, 1:3), an array
%   of triangular patches t(:, 1:3), and an array of outer normal vectors
%   normals(:, 1:3). 
%   
%   Mesh generator for planar shapes: Copyright DISTMESH 2004-2012 Per-Olof
%   Persson. 
%   GUI: Mr. Xingchi Dai for NEVA EM
%   Algorithm and its implementation: Sergey N Makarov
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2105, 1st ed.

%%  EM constants
eps0        = 8.85418782e-012;                  %   dielectric permittivity of vacuum(~air)
mu0         = 1.25663706e-006;                  %   magnetic permeability of vacuum(~air)
c0          = 1/sqrt(eps0*mu0);                 %   speed of light in vacuum(~air)
eta0        = sqrt(mu0/eps0);                   %   vacuum/air impedance, ~377 ohms

%%  GUI window - input parameters
global P;
global t;
global c;
global Area;
global Center;
global normals;
global Size;
global Points;
global ProjectFile;
ProjectFile = [];

%   Parameters for conducting objects (uimenu)
objecttype{1}           = 'circle';       %   'none' or 'plate', 'sphere', 'brick', or 'cylinder' is allowed

%   Data for the plate (if 'plate' is selected as a first object) (uitable)
strge.L(1)              = 0.2;         %   plate length, m
strge.W(1)              = 0.2;         %   plate width, m
strge.X(1)              = 0;           %   center position x, m
strge.Y(1)              = 0;           %   center position y, m
strge.Z(1)              = 0;           %   center position z, m
strge.ax(1)             = 0;           %   rotate about the x-axis, deg
strge.ay(1)             = 0;           %   rotate about the y-axis, deg
strge.az(1)             = 0;           %   rotate about the z-axis, deg
strge.par(1)            = 0.25;        %   uniform (0) or non-uniform (1) grid
strge.Tr(1)             = 300;

%   Data for the sphere (if 'sphere' is selected as a first object) (uitable)
strsphere.a(1)          = 0.04;        %   sphere radius, m
strsphere.X(1)          = 0;           %   center position x, m
strsphere.Y(1)          = 0;           %   center position y, m
strsphere.Z(1)          = 0;           %   center position z, m
strsphere.Tr(1)  = 300;                %   approximate number of triangular patches

%   Data for the brick (if 'brick' is selected as a first object) (uitable)
strbrick.L(1)          = 0.20;         %    brick length, m
strbrick.W(1)          = 0.20;         %    brick width, m
strbrick.H(1)          = 0.02;         %    brick height, m
strbrick.X(1)          = 0;            %    center position x, m
strbrick.Y(1)          = 0;            %    center position y, m
strbrick.Z(1)          = 0;            %    center position z, m
strbrick.ax(1)         = 0;            %    rotate about the x-axis, deg
strbrick.ay(1)         = 0;            %    rotate about the y-axis, deg
strbrick.az(1)         = 0;            %    rotate about the z-axis, deg
strbrick.par(1)  = 0.25;               %    uniform (0) or non-uniform (1) grid
strbrick.Tr(1)  = 300;

%   Data for the cylinder (if 'cylinder' is selected as a first object) (uitable)
strcylinder.R(1)          = 0.01;       %     cylinder radius, m
strcylinder.H(1)          = 0.05;       %     cylinder height, m
strcylinder.X(1)          = 0;          %     center position x, m
strcylinder.Y(1)          = 0;          %     center position y, m
strcylinder.Z(1)          = 0.;         %     center position z, m
strcylinder.ax(1)         = 90;         %     rotate about the x-axis, deg
strcylinder.ay(1)         = 0;          %     rotate about the y-axis, deg
strcylinder.az(1)         = 0;          %     rotate about the z-axis, deg
strcylinder.Tr(1)         = 6;
strcylinder.par(1)        = 1;

% Data for the circle (if 'circle' is selected as a first object) (uitable)
strcircle.R(1)              = 0.1;         %   circle radius, m
strcircle.X(1)              = 0;           %   center position x, m
strcircle.Y(1)              = 0;           %   center position y, m
strcircle.Z(1)              = 0.0;         %   center position z, m
strcircle.ax(1)             = 0;           %   rotate about the x-axis, deg
strcircle.ay(1)             = 0;           %   rotate about the y-axis, deg
strcircle.az(1)             = 0;           %   rotate about the z-axis, deg
strcircle.par(1)            = 0.25;        %   uniform (0) or non-uniform (1) grid
strcircle.Tr(1)             = 300;

strout.PatchesTotal = 0;
strout.NodesTotal = 0;
strout.quality = 0;
strout.averagequality = 0;
strout.time = 0;       

%%  End of GUI window - input parameters

%%  End of GUI window - output parameters

%   Input Parameters--------
%   Initialize handles structure
handles = struct();
handles.geometry = 0;   %   The real handles
handles.simulate = 0;   %   Only an indicator that the simulations are complete
handles.object   = 2;
%   Create all UI controls
build_gui();
geometry();

%%  Functions

    function cleaning
        %   cleaning figure without re-plotting
        io = 0;
        if handles.geometry ~= 0; delete(handles.geometry); handles.geometry = 0; end
        handles.geometry = subplot(1,1,1,'Parent',handles.uipanel2);
        graphics_S23(io, P, t, Center, normals);
    end

    function geometry()        
        %%   Conducting body mesh (it will be a cell array here)
        time1 = cputime;
        for m = 1 % size(objecttype, 2)
            if strcmp(objecttype{m}, 'plate')
                [PD{m}, tD{m}] = plate(strge.L(m), strge.W(m), strge.Tr(m), strge.par(m), 0, 0, 1);
                [PD{m}] = rotatex(PD{m}, strge.ax(m));
                [PD{m}] = rotatey(PD{m}, strge.ay(m));
                [PD{m}] = rotatez(PD{m}, strge.az(m));
                [centersD{m}, normalsD{m}] = normcenters(PD{m}, tD{m});
                PD{m}(:, 1) = PD{m}(:, 1) + strge.X(m);
                PD{m}(:, 2) = PD{m}(:, 2) + strge.Y(m);
                PD{m}(:, 3) = PD{m}(:, 3) + strge.Z(m);
                centersD{m}(:, 1) = centersD{m}(:, 1) + strge.X(m);
                centersD{m}(:, 2) = centersD{m}(:, 2) + strge.Y(m);
                centersD{m}(:, 3) = centersD{m}(:, 3) + strge.Z(m);               
            end
            if strcmp(objecttype{m}, 'sphere')
                [PD{m}, tD{m}] = sphere(1, strsphere.Tr(m));
                PD{m} = PD{m}*strsphere.a(m);
                [centersD{m}, normalsD{m}] = normcenters(PD{m}, tD{m});
                PD{m}(:, 1) = PD{m}(:, 1) + strsphere.X(m);
                PD{m}(:, 2) = PD{m}(:, 2) + strsphere.Y(m);
                PD{m}(:, 3) = PD{m}(:, 3) + strsphere.Z(m);            
            end
            if strcmp(objecttype{m}, 'brick')
                [PD{m}, tD{m}] = brick(strbrick.L(m), strbrick.W(m), strbrick.H(m), strbrick.Tr(m), strbrick.par(m));
                [PD{m}] = rotatex(PD{m}, strbrick.ax(m));
                [PD{m}] = rotatey(PD{m}, strbrick.ay(m));
                [PD{m}] = rotatez(PD{m}, strbrick.az(m));
                [centersD{m}, normalsD{m}] = normcenters(PD{m}, tD{m});
                PD{m}(:, 1) = PD{m}(:, 1) + strbrick.X(m);
                PD{m}(:, 2) = PD{m}(:, 2) + strbrick.Y(m);
                PD{m}(:, 3) = PD{m}(:, 3) + strbrick.Z(m);
                centersD{m}(:, 1) = centersD{m}(:, 1) + strbrick.X(m);
                centersD{m}(:, 2) = centersD{m}(:, 2) + strbrick.Y(m);
                centersD{m}(:, 3) = centersD{m}(:, 3) + strbrick.Z(m);              
            end
            if strcmp(objecttype{m}, 'cylinder')
                [PD{m}, tD{m}]  = cylinder(strcylinder.R(m), strcylinder.H(m), strcylinder.Tr(m), strcylinder.par(m));
                [PD{m}] = rotatex(PD{m}, strcylinder.ax(m));
                [PD{m}] = rotatey(PD{m}, strcylinder.ay(m));
                [PD{m}] = rotatez(PD{m}, strcylinder.az(m));
                [centersD{m}, normalsD{m}] = normcenters(PD{m}, tD{m});
                PD{m}(:, 1) = PD{m}(:, 1) + strcylinder.X(m);
                PD{m}(:, 2) = PD{m}(:, 2) + strcylinder.Y(m);
                PD{m}(:, 3) = PD{m}(:, 3) + strcylinder.Z(m);
                centersD{m}(:, 1) = centersD{m}(:, 1) + strcylinder.X(m);
                centersD{m}(:, 2) = centersD{m}(:, 2) + strcylinder.Y(m);
                centersD{m}(:, 3) = centersD{m}(:, 3) + strcylinder.Z(m);          
            end
            if strcmp(objecttype{m}, 'circle')
                [PD{m}, tD{m}] = circle(strcircle.R(m), strcircle.Tr(m),strcircle.par(m), 0, 0, 1);
                [PD{m}] = rotatex(PD{m}, strcircle.ax(m));
                [PD{m}] = rotatey(PD{m}, strcircle.ay(m));
                [PD{m}] = rotatez(PD{m}, strcircle.az(m));
                [centersD{m}, normalsD{m}] = normcenters(PD{m}, tD{m});
                PD{m}(:, 1) = PD{m}(:, 1) + strcircle.X(m);
                PD{m}(:, 2) = PD{m}(:, 2) + strcircle.Y(m);
                PD{m}(:, 3) = PD{m}(:, 3) + strcircle.Z(m);
                centersD{m}(:, 1) = centersD{m}(:, 1) + strcircle.X(m);
                centersD{m}(:, 2) = centersD{m}(:, 2) + strcircle.Y(m);
                centersD{m}(:, 3) = centersD{m}(:, 3) + strcircle.Z(m);
            end
        end
        time1 = cputime - time1;
        %%  Combining all meshes together (while keeping the fourth index)
        P = [];
        t = [];
        normals = [];
        for m = 1 
            tD{m}(:, 1:3)   = tD{m}(:, 1:3) + size(P, 1);
            P               = [P' PD{m}']';
            t               = [t' tD{m}']';
            normals         = [normals' normalsD{m}']';     %   normals
        end        
        Center = (P(t(:, 1),:)+P(t(:, 2),:)+P(t(:, 3),:))/3;
        %%  Input graphics
        strout.NodesTotal   = size(P, 1);  
        strout.PatchesTotal = size(t, 1);        
        strout.time = time1;               
        strout.quality          = min(simpqual(P, t));
        strout.averagequality   = mean(simpqual(P, t));
        set(0, 'CurrentFigure', handles.figure1);
        if handles.geometry ~= 0; delete(handles.geometry); handles.geometry = 0; end
        handles.geometry = subplot(1,1,1,'Parent',handles.uipanel2);
        io = 0;
        graphics_S23(io, P, t, Center, normals);
    end

    function simulate()             
        save ('mesh', 'P', 't', 'normals');
    end

%% ---------------------------------------------------------------------------
    function build_gui()
        % Creation of all uicontrols
        %---Basic Position parameters ---
        
        %% global Position Parameters
        %global Position_local;  Position_local  = [0.2, 0.27, 2.3/12, 6/12];
        global Position_b1;     Position_b1     = [1/12, 0.5/12, 4/12, 1/12];
        global Position_b2;     Position_b2     = [7/12, 0.5/12, 4/12, 1/12];
        global Position_b3;     Position_b3     = [0.1/12, 10.5/12,11.9/12, 1.5/12];
        global Position_table;  Position_table  = [0.1/12, 2/12, 1, 8/12];
        Position_units = 'pixels';
        
        % --- FIGURE -------------------------------------
        handles.figure1 = figure( ...
            'Tag', 'figure1', ...
            'Units', 'characters', ...
            'Position', [103.8 23.0769230769231 135.6 38.4615384615385], ...
            'Name', 'E00 - Mesh generator for simple shapes: uses DISTMESH 2004-2012 Per-Olof Persson', ...
            'MenuBar', 'none', ...
            'NumberTitle', 'off', ...
            'Color', get(0,'DefaultUicontrolBackgroundColor'), ...
            'Resize', 'on',...
            'Toolbar','figure');
        set( handles.figure1, ...
            'Units', 'pixels' );
        
        %get your display size
        screenSize = get(0, 'ScreenSize');
        
        %calculate the center of the display
        position = get( handles.figure1, ...
            'Position' );
        position(2) = (screenSize(4)-position(4))/2;
        position(1) = (screenSize(3)-position(3))/2;
        %center the window
        set( handles.figure1, ...
            'Position', position );
        
        
        Position_local = get(handles.figure1,'position');
        Position_local(1) = 0.5*Position_local(1) ;
        Position_local(2) = 0.87*Position_local(2) ;
        Position_local(3) = 0.48*Position_local(3);
        Position_local(4) = Position_local(4);
        
        
        hToolLegend1 = findall(handles.figure1,'tag','Standard.NewFigure');
        hToolLegend2 = findall(handles.figure1,'tag','Standard.FileOpen');
        hToolLegend3 = findall(handles.figure1,'tag','Standard.SaveFigure');
        hToolLegend4 = findall(handles.figure1,'tag','Standard.EditPlot');
        hToolLegend5 = findall(handles.figure1,'tag','Exploration.Pan');
        hToolLegend6 = findall(handles.figure1,'tag','DataManager.Linking');
        hToolLegend7 = findall(handles.figure1,'tag','Annotation.InsertLegend');
        hToolLegend8 = findall(handles.figure1,'tag','Exploration.Brushing');
        hToolLegend9 = findall(handles.figure1,'tag','Standard.PrintFigure');
        hToolLegend10 = findall(handles.figure1,'tag','Exploration.ZoomIn');
        hToolLegend11 = findall(handles.figure1,'tag','Exploration.ZoomOut');
        
        delete(hToolLegend1);
        delete(hToolLegend2);
        delete(hToolLegend4);
        delete(hToolLegend5);
        delete(hToolLegend6);
        delete(hToolLegend7);
        delete(hToolLegend8);
        delete(hToolLegend10);
        delete(hToolLegend11);
        
        
        % Make a Print-view
        hToolbar = findall(gcf,'tag','FigureToolBar');
        hPrintButton = findall(hToolbar,'tag','Standard.PrintFigure');
        set(hPrintButton, 'ClickedCallback','printpreview(gcbf)', ...
            'TooltipString','Print Preview');
        
        
        % --- Menus --------------------------------------
        
        handles.fFileMenu      =   uimenu(...       % File menu
            'Parent', handles.figure1,...
            'HandleVisibility','callback', ...
            'Label','Project');
        
        uimenu(...
            'Parent', handles.fFileMenu,...
            'HandleVisibility','callback', ...
            'Label','New Project',...
            'Callback', {@New_Callback});
        
        
        function New_Callback(~, ~)
            ProjectFile = 'project1';
        end
        
        uimenu(...
            'Parent', handles.fFileMenu,...
            'HandleVisibility','callback', ...
            'Label','Open',...
            'Callback', {@Open_Callback});
        
        function Open_Callback(~, ~)
            
            file = uigetfile('*.mat', 'Select project file to open');            
            if ~isequal(file, 0) 
                ProjectFile = file;  
                M   = load(file);             
                strge               = M.strge;
                strsphere           = M.strsphere;
                strbrick            = M.strbrick;
                strcylinder         = M.strcylinder;
                strcircle           = M.strcircle;
                objecttype          = M.objecttype;
                switch objecttype{1}
                    case 'plate'
                        handles.object = 1;
                    case 'circle'
                        handles.object = 2;
                    case 'sphere'
                        handles.object = 3;
                    case 'brick'
                        handles.object = 4;
                    case 'cylinder'
                        handles.object = 5;
                end                
                geometry();
                        
            end
        end
        
        uimenu(...
            'Parent', handles.fFileMenu,...
            'HandleVisibility','callback', ...
            'Label','Save',...
            'Callback', {@Save_Callback});
        
        function Save_Callback(~, ~)
            if ~isempty(ProjectFile)
                save(ProjectFile);
            else
                save('E00project');
            end
        end
         
        uimenu(...
            'Parent', handles.fFileMenu,...
            'HandleVisibility','callback', ...
            'Label','Save As',...
            'Callback', {@SaveAs_Callback});
        
        function SaveAs_Callback(~,~)
            file = uiputfile('*.mat','Save project file as');
            if ~isequal(file, 0)
                save(file);
            end
        end
            
        handles.fFigureMenu      =   uimenu(...        % File menu
            'Parent', handles.figure1,...
            'HandleVisibility','callback', ...
            'Label','Save figure',...
            'callback',@SaveFigureCallback);
        
        function SaveFigureCallback(~,~)
            file = uiputfile({'*.fig';'*.png';'*.jpeg'},'Save figure as');
            if ~isequal(file, 0)
                if handles.geometry ~= 0; saveas(handles.geometry,file);end
                if handles.simulate ~= 0; saveas(handles.simulate,file);end
            end
        end
                        
        handles.fFirstMenu   =   uimenu(...             % File menu
            'Parent', handles.figure1,...
            'HandleVisibility','callback', ...
            'Label','Create mesh');
        
        handles.fObjectTypeMenu      =   uimenu(...     % File menu
            'Parent', handles.fFirstMenu,...
            'HandleVisibility','callback', ...
            'Label','Object Type',...
            'Callback', @fObjectTypeCallback);
        
        
        function fObjectTypeCallback(~,~)
            % --- RADIO BUTTONS -------------------------------------
            handles.figure2 = figure( ...
                'Tag', 'figure1', ...
                'Units', 'characters', ...
                'Position', [103.8 38.6153846153846 28.4 15.7692307692308], ...
                'Name', 'Object Type', ...
                'MenuBar', 'none', ...
                'NumberTitle', 'off', ...
                'Color', get(0,'DefaultUicontrolBackgroundColor'), ...
                'Resize', 'on');
            
            handles.uipanel3 = uibuttongroup( ...
                'Parent', handles.figure2, ...
                'Tag', 'uipanel1', ...
                'Units', 'normalized', ...
                'Position', [0.0141823161189358 0.0334928229665072 0.965277777777778 0.961722488038277], ...
                'FontSize', 10.6666666666667, ...
                'FontUnits', 'pixels', ...
                'Title', 'Object Type');
            
            set(handles.uipanel3,'SelectionChangeFcn',@uibuttongroup1_SelectionChangeFcn);
            
            
            function uibuttongroup1_SelectionChangeFcn(~,eventdata)
                switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
                    case 'plate'
                        handles.object = 1;
                        objecttype{1} = 'plate';
                        set(0, 'CurrentFigure', handles.figure1);
                        geometry();              
                    case 'circle'
                        handles.object = 2;
                        objecttype{1} = 'circle';
                        set(0, 'CurrentFigure', handles.figure1);
                        geometry();               
                    case 'sphere'
                        handles.object = 3;
                        objecttype{1} = 'sphere';
                        set(0, 'CurrentFigure', handles.figure1);
                        geometry();           
                    case 'brick'
                        handles.object = 4;
                        objecttype{1} = 'brick';
                        set(0, 'CurrentFigure', handles.figure1);
                        geometry();          
                    case 'cylinder'
                        handles.object = 5;
                        objecttype{1} = 'cylinder';
                        set(0, 'CurrentFigure', handles.figure1);
                        geometry();
                    otherwise
                        errordlg('Fatal error! Contact us!');
                end
            end
                        
            handles.radiobutton7 = uicontrol( ...
                'Parent', handles.uipanel3, ...
                'Tag', 'plate', ...
                'Style', 'radiobutton', ...
                'Units', 'normalized', ...
                'Position', [0.165413533834586 0.877777777777778 0.654135338345865 0.127777777777778], ...
                'FontSize', 10.6666666666667, ...
                'FontUnits', 'pixels', ...
                'String', 'Plate');
            
            
            handles.radiobutton8 = uicontrol( ...
                'Parent', handles.uipanel3, ...
                'Tag', 'circle', ...
                'Style', 'radiobutton', ...
                'Units', 'normalized', ...
                'Position', [0.165413533834586 0.666666666666667 0.654135338345865 0.127777777777778], ...
                'FontSize', 10.6666666666667, ...
                'FontUnits', 'pixels', ...
                'String', 'Circle');
                        
            
            handles.radiobutton9 = uicontrol( ...
                'Parent', handles.uipanel3, ...
                'Tag', 'sphere', ...
                'Style', 'radiobutton', ...
                'Units', 'normalized', ...
                'Position', [0.165413533834586 0.455555555555556 0.654135338345865 0.127777777777778], ...
                'FontSize', 10.6666666666667, ...
                'FontUnits', 'pixels', ...
                'String', 'Sphere');
            
            
            handles.radiobutton10 = uicontrol( ...
                'Parent', handles.uipanel3, ...
                'Tag', 'brick', ...
                'Style', 'radiobutton', ...
                'Units', 'normalized', ...
                'Position', [0.165413533834586 0.244444444444444 0.654135338345865 0.127777777777778], ...
                'FontSize', 10.6666666666667, ...
                'FontUnits', 'pixels', ...
                'String', 'Brick');
            
            
            handles.radiobutton11 = uicontrol( ...
                'Parent', handles.uipanel3, ...
                'Tag', 'cylinder', ...
                'Style', 'radiobutton', ...
                'Units', 'normalized', ...
                'Position', [0.165413533834586 0.0333333333333333 0.654135338345865 0.127777777777778], ...
                'FontSize', 10.6666666666667, ...
                'FontUnits', 'pixels', ...
                'String', 'Cylinder');
            
            set(handles.uipanel3,'Selectedobject',handles.radiobutton10);
            
            switch handles.object
                case 1
                    set(handles.uipanel3,'Selectedobject',handles.radiobutton7);
                case 2
                    set(handles.uipanel3,'Selectedobject',handles.radiobutton8);
                case 3
                    set(handles.uipanel3,'Selectedobject',handles.radiobutton9);
                case 4
                    set(handles.uipanel3,'Selectedobject',handles.radiobutton10);
                case 5
                    set(handles.uipanel3,'Selectedobject',handles.radiobutton11);
                otherwise
                    errordlg('Fatal error;Contact your Administrator');
            end
            
        end
                
        handles.fParametersMenu      =   uimenu(...       % File menu
            'Parent', handles.fFirstMenu,...
            'HandleVisibility','callback', ...
            'Label','Parameters',...
            'Callback', @fPreParametersCallback);
        
        function fPreParametersCallback(~,~)
            
            switch handles.object
                case 1
                    fConductingPlateCallback;
                case 2
                    fConductingCircleCallback;
                case 3
                    fConductingSphereCallback;
                case 4
                    fConductingBrickCallback;
                case 5
                    fConductingCylinderCallback;
                otherwise
                    errordlg('Fatal error! Contact us!');
            end            
            
        end
        
        
        function fConductingPlateCallback(~, ~)
            strge0 = strge;
            % Callback function run when_______________________________
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            local = figure('Units', Position_units, 'Position', Position_local);
            set(local, 'Name', 'Conductor setting', 'NumberTitle', 'off', 'MenuBar', 'none');
            b1 = uicontrol('Style', 'Pushbutton', 'Units', 'normalized', 'String', 'Update','Position', Position_b1, 'Parent', local, 'Callback', @B1_Callback);
            b2 = uicontrol('Style', 'Pushbutton', 'Units', 'normalized', 'String', 'Cancel', 'Position', Position_b2, 'Parent', local, 'Callback', @B2_Callback);
            b3 = uicontrol('Style', 'Text', 'Units', 'normalized', 'String', {'parameters of the plate setup';'Press ENTER after changing each parameter'}, 'Position', Position_b3, 'Parent', local, 'BackgroundColor', [0.8, 0.8, 0.8]);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            data = {...
                'Plate length,m',strge.L(1);...
                'Plate width,m',strge.W(1);...
                'Center position x,m',strge.X(1);...
                'Center position y,m',strge.Y(1);...
                'Center position z,m',strge.Z(1);...
                'Rotate about x-axis, deg',strge.ax(1);...
                'Rotate about y-axis, deg',strge.ay(1);...
                'Rotate about z-axis, deg',strge.az(1);...
                'Approximate number of triangles',strge.Tr(1);...
                'Triangle size at boundaries vs. center size', strge.par(1)};
            ColumnName =   {'Name',   'Value'};
            columnformat    = {'char'};
            columneditable  =  [true];
            uitable('Units', 'normalized', 'Position', Position_table,...
                'Data', data, 'Visible', 'on', 'BackgroundColor', [1 1 1],...
                'RowName', [], 'ColumnName', ColumnName, 'ColumnFormat', columnformat, ...
                'ColumnEditable', columneditable, 'ToolTipString', 'You have chosen the plate case',...
                'ColumnWidth',{200,100},...
                'CellEditCallback', {@edit_callback});
            function edit_callback(hObject, ~)
                data = get(hObject,'Data');
                strge.L(1)         = data{11};
                strge.W(1)         = data{12};
                strge.X(1)         = data{13};
                strge.Y(1)         = data{14};
                strge.Z(1)         = data{15};
                strge.ax(1)        = data{16};
                strge.ay(1)        = data{17};
                strge.az(1)        = data{18};
                strge.Tr(1)        = data{19};
                strge.par(1)       = data{20};
            end
            
            function B1_Callback(~, ~, ~)
                set(0, 'CurrentFigure', handles.figure1);
                geometry();           
                delete(local);
            end
            function B2_Callback(~, ~, ~)
                strge = strge0;
                delete(local);
            end
            
        end
        
        %% Circle case
        function fConductingCircleCallback(~, ~)
            strcircle0 = strcircle;
            % Callback function run when_______________________________
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            local = figure('Units', Position_units, 'Position', Position_local);
            set(local, 'Name', 'Conductor setting', 'NumberTitle', 'off', 'MenuBar', 'none');
            b1 = uicontrol('Style', 'Pushbutton', 'Units', 'normalized', 'String', 'Update','Position', Position_b1, 'Parent', local, 'Callback', @B1_Callback);
            b2 = uicontrol('Style', 'Pushbutton', 'Units', 'normalized', 'String', 'Cancel', 'Position', Position_b2, 'Parent', local, 'Callback', @B2_Callback);
            b3 = uicontrol('Style', 'Text', 'Units', 'normalized', 'String', {'parameters of circle setup';'Press ENTER after changing each parameter'}, 'Position', Position_b3, 'Parent', local, 'BackgroundColor', [0.8, 0.8, 0.8]);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            data = {...
                'Circle radius,m',strcircle.R(1);...
                'Center position x,m',strcircle.X(1);...
                'Center position y,m',strcircle.Y(1);...
                'Center position z,m',strcircle.Z(1);...
                'Rotate about x-axis, deg',strcircle.ax(1);...
                'Rotate about y-axis, deg',strcircle.ay(1);...
                'Rotate about z-axis, deg',strcircle.az(1);...
                'Approximate numbers of triangles',strcircle.Tr(1);...
                'Triangle size at boundaries vs. center size', strcircle.par(1)};
            ColumnName =   {'Name',   'Value'};
            columnformat    = {'char'};
            columneditable  =  [true];
            uitable('Units', 'normalized', 'Position', Position_table,...
                'Data', data, 'Visible', 'on', 'BackgroundColor', [1 1 1],...
                'RowName', [], 'ColumnName', ColumnName, 'ColumnFormat', columnformat, ...
                'ColumnEditable', columneditable, 'ToolTipString', 'You have chosen the circle',...
                'ColumnWidth',{200,100},...
                'CellEditCallback', {@edit_callback});
            function edit_callback(hObject, ~)
                data = get(hObject,'Data');
                strcircle.R(1)         = data{10};
                strcircle.X(1)         = data{11};
                strcircle.Y(1)         = data{12};
                strcircle.Z(1)         = data{13};
                strcircle.ax(1)        = data{14};
                strcircle.ay(1)        = data{15};
                strcircle.az(1)        = data{16};
                strcircle.Tr(1)         = data{17};
                strcircle.par(1) = data{18};
                                
            end
            
            function B1_Callback(~, ~, ~)
                set(0, 'CurrentFigure', handles.figure1);
                geometry();
                delete(local);
            end
            function B2_Callback(~, ~, ~)
                strcircle = strcircle0;
                delete(local);
            end
            
        end
        %% Sphere case
        
        function fConductingSphereCallback(~, ~)
            strsphere0 = strsphere;
            % Callback function run when_______________________________
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            local = figure('Units', Position_units, 'Position', Position_local);
            set(local, 'Name', 'Conductor setting', 'NumberTitle', 'off', 'MenuBar', 'none');
            b1 = uicontrol('Style', 'Pushbutton', 'Units', 'normalized', 'String', 'Update','Position', Position_b1, 'Parent', local, 'Callback', @B1_Callback);
            b2 = uicontrol('Style', 'Pushbutton', 'Units', 'normalized', 'String', 'Cancel', 'Position', Position_b2, 'Parent', local, 'Callback', @B2_Callback);
            b3 = uicontrol('Style', 'Text', 'Units', 'normalized', 'String', {'parameters of sphere setup';'Press ENTER after changing each parameter'}, 'Position', Position_b3, 'Parent', local, 'BackgroundColor', [0.8, 0.8, 0.8]);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            data = {...
                'Sphere radius, m',strsphere.a(1);...
                'Center position x,m',strsphere.X(1);...
                'Center position y,m',strsphere.Y(1);...
                'Center position z,m',strsphere.Z(1);...
                'Approx.# of triangular patches',strsphere.Tr(1)};
            ColumnName =   {'Name',   'Value'};
            columnformat    = {'char'};
            columneditable  =  [true];
            uitable('Units', 'normalized', 'Position', Position_table,...
                'Data', data, 'Visible', 'on', 'BackgroundColor', [1 1 1],...
                'RowName', [], 'ColumnName', ColumnName, 'ColumnFormat', columnformat, ...
                'ColumnEditable', columneditable, 'ToolTipString', 'You have chosen the sphere case ',...
                'ColumnWidth',{200,100},...
                'CellEditCallback', {@edit_callback});
            function edit_callback(hObject, ~)
                data = get(hObject,'Data');
                strsphere.a(1)                 = data{6};
                strsphere.X(1)                 = data{7};
                strsphere.Y(1)                 = data{8};
                strsphere.Z(1)                 = data{9};
                strsphere.Tr(1)                = data{10};
            end
            
            function B1_Callback(~, ~, ~)
                set(0, 'CurrentFigure', handles.figure1);
                geometry();
                delete(local);
            end
            function B2_Callback(~, ~, ~)
                strsphere = strsphere0;
                delete(local);
            end
            
        end
        %% Brick Case
        
        function fConductingBrickCallback(~, ~)
            strbrick0 = strbrick;
            % Callback function run when_______________________________
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            local = figure('Units', Position_units, 'Position', Position_local);
            set(local, 'Name', 'Conductor setting', 'NumberTitle', 'off', 'MenuBar', 'none');
            b1 = uicontrol('Style', 'Pushbutton', 'Units', 'normalized', 'String', 'Update','Position', Position_b1, 'Parent', local, 'Callback', @B1_Callback);
            b2 = uicontrol('Style', 'Pushbutton', 'Units', 'normalized', 'String', 'Cancel', 'Position', Position_b2, 'Parent', local, 'Callback', @B2_Callback);
            b3 = uicontrol('Style', 'Text', 'Units', 'normalized', 'String', {'Parameters of brick setup';'Press ENTER after changing each parameter'}, 'Position', Position_b3, 'Parent', local, 'BackgroundColor', [0.8, 0.8, 0.8]);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            data = {...
                'Brick length,m',strbrick.L(1);...
                'Brick width,m',strbrick.W(1);...
                'Brick height,m',strbrick.H(1);...
                'Center position x,m',strbrick.X(1);...
                'Center position y,m',strbrick.Y(1);...
                'Center position z,m',strbrick.Z(1);...
                'Rotate about x-axis, deg',strbrick.ax(1);...
                'Rotate about y-axis, deg',strbrick.ay(1);...
                'Rotate about z-axis, deg',strbrick.az(1);...
                'Approximate # of triangles - top face',strbrick.Tr(1);...
                'Triangle size at boundaries vs. center size', strbrick.par(1)};
            ColumnName =   {'Name',   'Value'};
            columnformat    = {'char'};
            columneditable  =  [true];
            uitable('Units', 'normalized', 'Position', Position_table,...
                'Data', data, 'Visible', 'on', 'BackgroundColor', [1 1 1],...
                'RowName', [], 'ColumnName', ColumnName, 'ColumnFormat', columnformat, ...
                'ColumnEditable', columneditable, 'ToolTipString', 'You have chosen the brick case',...
                'ColumnWidth',{200,100},...
                'CellEditCallback', {@edit_callback});
            function edit_callback(hObject, ~)
                data = get(hObject,'Data');
                strbrick.L(1)         = data{12};
                strbrick.W(1)         = data{13};
                strbrick.H(1)         = data{14};
                strbrick.X(1)         = data{15};
                strbrick.Y(1)         = data{16};
                strbrick.Z(1)         = data{17};
                strbrick.ax(1)        = data{18};
                strbrick.ay(1)        = data{19};
                strbrick.az(1)        = data{20};
                strbrick.Tr(1)        = data{21};
                strbrick.par(1)       = data{22};
            end
            
            function B1_Callback(~, ~, ~)
                set(0, 'CurrentFigure', handles.figure1);
                geometry();
                delete(local);
            end
            function B2_Callback(~, ~, ~)
                strbrick = strbrick0;
                delete(local);
            end
            
        end
        
        %% Cylinder case
        
        function fConductingCylinderCallback(~, ~)
            strcylinder0 = strcylinder;
            % Callback function run when_______________________________
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            local = figure('Units', Position_units, 'Position', Position_local);
            set(local, 'Name', 'Conductor setting', 'NumberTitle', 'off', 'MenuBar', 'none');
            b1 = uicontrol('Style', 'Pushbutton', 'Units', 'normalized', 'String', 'Update','Position', Position_b1, 'Parent', local, 'Callback', @B1_Callback);
            b2 = uicontrol('Style', 'Pushbutton', 'Units', 'normalized', 'String', 'Cancel', 'Position', Position_b2, 'Parent', local, 'Callback', @B2_Callback);
            b3 = uicontrol('Style', 'Text', 'Units', 'normalized', 'String', {'parameters of cylinder setup';'Press ENTER after changing each parameter'}, 'Position', Position_b3, 'Parent', local, 'BackgroundColor', [0.8, 0.8, 0.8]);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            data = {...
                'Cylinder radius,m',strcylinder.R(1);...
                'Cylinder height,m',strcylinder.H(1);...
                'Center position x,m',strcylinder.X(1);...
                'Center position y,m',strcylinder.Y(1);...
                'Center position z,m',strcylinder.Z(1);...
                'Rotate about x-axis, deg',strcylinder.ax(1);...
                'Rotate about y-axis, deg',strcylinder.ay(1);...
                'Rotate about z-axis, deg',strcylinder.az(1);...
                'Approximate # of triangles - top face',strcylinder.Tr(1);...
                'Triangle size at boundaries vs. center size',strcylinder.par(1)};
            ColumnName =   {'Name',   'Value'};
            columnformat    = {'char'};
            columneditable  =  [true];
            uitable('Units', 'normalized', 'Position', Position_table,...
                'Data', data, 'Visible', 'on', 'BackgroundColor', [1 1 1],...
                'RowName', [], 'ColumnName', ColumnName, 'ColumnFormat', columnformat, ...
                'ColumnEditable', columneditable, 'ToolTipString', 'You have chosen the cylinder',...
                'ColumnWidth',{200,100},...
                'CellEditCallback', {@edit_callback});
            function edit_callback(hObject, ~)
                data = get(hObject,'Data');
                strcylinder.R(1)         = data{11};
                strcylinder.H(1)         = data{12};
                strcylinder.X(1)         = data{13};
                strcylinder.Y(1)         = data{14};
                strcylinder.Z(1)         = data{15};
                strcylinder.ax(1)        = data{16};
                strcylinder.ay(1)        = data{17};
                strcylinder.az(1)        = data{18};
                strcylinder.Tr(1)        = data{19};
                strcylinder.par(1)       = data{20};
            end
            
            function B1_Callback(~, ~, ~)
                set(0, 'CurrentFigure', handles.figure1);
                geometry();
                delete(local);
            end
            function B2_Callback(~, ~, ~)
                strcylinder = strcylinder0;
                delete(local);
            end
            
        end
        
        % --- Normals Menu ----
        handles.fNormalMenu    =   uimenu(...       % File menu
            'Parent',handles.figure1,...
            'HandleVisibility','callback', ...
            'Label','Show normal vectors',...
            'callback',@fNormalSetupCallback);
        
        function fNormalSetupCallback(~, ~)
            handles.simulate = 0;
             graphics_S23(1, P, t, Center, normals);
        end
        
         % --- Output Result ---
        handles.fOutputDataMenu    =   uimenu(...       % File menu
            'Parent',handles.figure1,...
            'HandleVisibility','callback', ...
            'Label','Mesh data',...
            'callback',@fOutputDataCallback);
        
        function fOutputDataCallback(~, ~)
            % Callback function run when_______________________________
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            local = figure('Units', Position_units, 'Position', Position_local);
            set(local, 'Name', 'Output Data', 'NumberTitle', 'off', 'MenuBar', 'none');
            b2 = uicontrol('Style', 'Pushbutton', 'Units', 'normalized', 'String', 'Exit', 'Position', Position_b2, 'Parent', local, 'Callback', @B2_Callback);
            b3 = uicontrol('Style', 'Text', 'Units', 'normalized', 'String', 'Mesh data', 'Position', Position_b3, 'Parent', local, 'BackgroundColor', [0.8, 0.8, 0.8]);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ColumnName =   {'Name',   'Value'};
            columneditable  =  [true];
            data = {...
                'Total number of triangular patches',strout.PatchesTotal;...
                'Total number of nodes',              strout.NodesTotal;...
                'Minimum triangle quality in the mesh', strout.quality;...
                'Average triangle quality in the mesh', strout.averagequality;...
                'CPU time in sec for creating the mesh, s', strout.time};                        
            
            handles.outputtable = uitable('Units', 'normalized', 'Position', Position_table,...
                'Data',data,'Visible', 'on', 'BackgroundColor', [1 1 1],...
                'RowName', [], 'ColumnName', ColumnName, ...
                'ColumnEditable', columneditable, 'ToolTipString', 'Mesh data',...
                'ColumnWidth',{200,100});
            function B2_Callback(~, ~, ~)
                delete(local);
            end
        end
   
     % --- View Direction ---
        handles.fViewMenu    =   uimenu(...       % File menu
            'Parent',handles.figure1,...
            'HandleVisibility','callback', ...
            'Label','View',...
            'callback',@fViewCallback);
        
        function fViewCallback(~,~)
            
        % --- RADIO BUTTONS -------------------------------------
             handles.figure3 = figure( ...
                'Tag', 'figure1', ...
                'Units', 'characters', ...
                'Position', [103.8 38.6153846153846 28.4 20], ...
                'Name', 'View direction', ...
                'MenuBar', 'none', ...
                'NumberTitle', 'off', ...
                'Color', get(0,'DefaultUicontrolBackgroundColor'), ...
                'Resize', 'on');
            set( handles.figure3, ...
                'Units', 'pixels' );
            
            %calculate the center of the display
            position = get(handles.figure1,'position');
            position(1) = position(1)*0.75;
            position(2) = position(2)*1.3;
            position(3) = position(3)*0.25;
            position(4) = position(4)*0.5;
            %center the window
            set( handles.figure3, ...
                'Position', position ); 
            
            handles.uipanel3 = uibuttongroup( ...
                'Parent', handles.figure3, ...
                'Tag', 'uipanel1', ...
                'Units', 'normalized', ...
                'Position', [0.0141823161189358 0.0334928229665072 0.965277777777778 0.961722488038277], ...
                'FontSize', 10.6666666666667, ...
                'FontUnits', 'pixels', ...
                'Title', 'View direction');
            
            set(handles.uipanel3,'SelectionChangeFcn',@uibuttongroup1_SelectionChangeFcn);
            
            
            function uibuttongroup1_SelectionChangeFcn(~,eventdata)
                switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
                    case 'Current'
                        set(0, 'CurrentFigure', handles.figure1);                 
                        view(-47, 20);
                    case 'Front'
                        set(0, 'CurrentFigure', handles.figure1);
                        view(90,0);
                    case 'Rear'
                        set(0, 'CurrentFigure', handles.figure1);
                        view(270,0);
                    case 'Top'
                        set(0, 'CurrentFigure', handles.figure1);
                        view(90,90);
                    case 'Bottom'
                        set(0, 'CurrentFigure', handles.figure1);
                        view(90,-90);
                    case 'Left'
                        set(0, 'CurrentFigure', handles.figure1);
                        view(0,0);
                    case 'Right'
                        set(0, 'CurrentFigure', handles.figure1);
                        view(180,0);
                    otherwise
                        errordlg('Fatal error! Contact us!');
                end
            end
            
            
            handles.radiobutton7 = uicontrol( ...
                'Parent', handles.uipanel3, ...
                'Tag', 'Default', ...
                'Style', 'radiobutton', ...
                'Units', 'normalized', ...
                'Position', [0.165413533834586 0.86 0.654135338345865 0.127777777777778], ...
                'FontSize', 10.6666666666667, ...
                'FontUnits', 'pixels', ...
                'String', 'Default');
            
            
            handles.radiobutton8 = uicontrol( ...
                'Parent', handles.uipanel3, ...
                'Tag', 'Front', ...
                'Style', 'radiobutton', ...
                'Units', 'normalized', ...
                'Position', [0.165413533834586 0.72 0.654135338345865 0.127777777777778], ...
                'FontSize', 10.6666666666667, ...
                'FontUnits', 'pixels', ...
                'String', 'Front');
                        
            handles.radiobutton9 = uicontrol( ...
                'Parent', handles.uipanel3, ...
                'Tag', 'Rear', ...
                'Style', 'radiobutton', ...
                'Units', 'normalized', ...
                'Position', [0.165413533834586 0.58 0.654135338345865 0.127777777777778], ...
                'FontSize', 10.6666666666667, ...
                'FontUnits', 'pixels', ...
                'String', 'Rear');
                        
            handles.radiobutton10 = uicontrol( ...
                'Parent', handles.uipanel3, ...
                'Tag', 'Top', ...
                'Style', 'radiobutton', ...
                'Units', 'normalized', ...
                'Position', [0.165413533834586 0.44 0.654135338345865 0.127777777777778], ...
                'FontSize', 10.6666666666667, ...
                'FontUnits', 'pixels', ...
                'String', 'Top');
            
            handles.radiobutton11 = uicontrol( ...
                'Parent', handles.uipanel3, ...
                'Tag', 'Bottom', ...
                'Style', 'radiobutton', ...
                'Units', 'normalized', ...
                'Position', [0.165413533834586 0.28 0.654135338345865 0.127777777777778], ...
                'FontSize', 10.6666666666667, ...
                'FontUnits', 'pixels', ...
                'String', 'Bottom');
            
            handles.radiobutton12 = uicontrol( ...
                'Parent', handles.uipanel3, ...
                'Tag', 'Left', ...
                'Style', 'radiobutton', ...
                'Units', 'normalized', ...
                'Position', [0.165413533834586 0.14 0.654135338345865 0.127777777777778], ...
                'FontSize', 10.6666666666667, ...
                'FontUnits', 'pixels', ...
                'String', 'Left');
                
            handles.radiobutton13 = uicontrol( ...
                'Parent', handles.uipanel3, ...
                'Tag', 'Right', ...
                'Style', 'radiobutton', ...
                'Units', 'normalized', ...
                'Position', [0.165413533834586 0.00 0.654135338345865 0.127777777777778], ...
                'FontSize', 10.6666666666667, ...
                'FontUnits', 'pixels', ...
                'String', 'Right');
        end
        
        % --- PANELS -------------------------------------
        handles.uipanel1 = uipanel( ...
            'Parent', handles.figure1, ...
            'Tag', 'uipanel1', ...
            'Units', 'normalized', ...
            'Position', [0.0162241887905605 0.822 0.967551622418879 0.156], ...
            'FontSize', 7.6666666666666, ...
            'FontUnits', 'pixels', ...
            'Title', 'Abstract');
        
        handles.uipanel2 = uipanel( ...
            'Parent', handles.figure1, ...
            'Tag', 'uipanel2', ...
            'Units', 'normalized', ...
            'Position', [0.0162241887905605 0.104 0.967551622418879 0.718], ...
            'FontSize', 7.6666666666666, ...
            'FontUnits', 'pixels', ...
            'Title', 'Figure');
        
        % --- ABSTRACT PANEL----------------------------------
        handles.abstract = uicontrol( ...
            'Parent', handles.uipanel1, ...
            'Tag', 'Abstract Description', ...
            'HorizontalAlignment', 'left',...
            'Style', 'text', ...
            'Units', 'normalized', ...
            'Position', [0 0.2 1 0.7], ...
            'FontSize', 8, ...
            'FontUnits', 'pixels', ...           
            'String', 'This module outputs uniform and non-uniform triangular meshes for a plate, circle, brick, sphere, or a cylinder into *.mat file mesh.mat. The mesh includes an array of nodes P(:, 1:3), array of triangular patches t(:, 1:3), and array of outer normal vectors normals(:, 1:3). Mesh generator: Copyright DISTMESH 2004-2012 Per-Olof Persson. GUI: Xingchi Dai for NEVA.');
        % --- PUSHBUTTONS -------------------------------------
        handles.pushbutton2 = uicontrol( ...
            'Parent', handles.figure1, ...
            'Tag', 'pushbutton2', ...
            'Style', 'pushbutton', ...
            'Units', 'normalized', ...
            'Position', [0.55 0.0271317829457364 0.105413105413105 0.0542635658914729], ...
            'FontSize', 7.6666666666666, ...
            'FontUnits', 'pixels', ...
            'String', 'Reload',...
            'callback',@GeometryCallback);
        
        function GeometryCallback(~,~)
            handles.simulate = 0;
            geometry();
            cleaning;
        end
        
        handles.pushbutton3 = uicontrol( ...
            'Parent', handles.figure1, ...
            'Tag', 'pushbutton3', ...
            'Style', 'pushbutton', ...
            'Units', 'normalized', ...
            'Position', [0.7 0.0271317829457364 0.106837606837607 0.0542635658914729], ...
            'FontSize', 7.6666666666666, ...
            'FontUnits', 'pixels', ...
            'String', 'Exit',...
            'Callback',@ExitCallback);
        
        function ExitCallback(~,~)
            
            button = questdlg('Save project before closing?');
            
            switch button
                case 'Yes',
                    if ~isempty(ProjectFile)
                        save(ProjectFile);
                    else                        
                        save('E00project');
                    end
                    close all;
                case'No',
                    close all;
                case 'Cancel',
                    quit cancel;
            end
            
        end
        
        handles.pushbutton1 = uicontrol( ...
            'Parent', handles.figure1, ...
            'Tag', 'pushbutton1', ...
            'Style', 'pushbutton', ...
            'Units', 'normalized', ...
            'Position', [0.4 0.0271317829457364 0.105413105413105 0.0542635658914729], ...
            'FontSize', 7.6666666666666, ...
            'FontUnits', 'pixels', ...
            'String', 'Save mesh',...
            'callback',@SimulateCallback);
        
        function SimulateCallback(~,~)
            simulate();
        end
        
        handles.pushbutton4 = uicontrol( ...
            'Parent', handles.figure1, ...
            'Tag', 'pushbutton3', ...
            'Style', 'pushbutton', ...
            'Units', 'normalized', ...
            'Position', [0.2 0.0271317829457364 0.15 0.0542635658914729], ...
            'FontSize', 7.6666666666666, ...
            'FontUnits', 'pixels', ...
            'String', 'Separate window',...
            'Callback',@OpenZoomCallback);
        
        function OpenZoomCallback(~,~)
            
            handles.figure2 = figure( ...
                'Tag', 'figure2', ...
                'Units', 'characters', ...
                'Position', [103.8 23.0769230769231 50 20], ...
                'Name', 'Figure window', ...
                'MenuBar', 'none', ...
                'NumberTitle', 'off', ...
                'Color', get(0,'DefaultUicontrolBackgroundColor'), ...
                'Resize', 'on',...
                'Toolbar','figure');
            
            set( handles.figure2, ...
                'Units', 'pixels' );
            
            %calculate the center of the display
            position = get(handles.figure1,'position');
            position(1) = position(1)*0.75;
            position(2) = position(2)*0.65;
            position(3) = position(3)*0.75;
            position(4) = position(4)*0.75;
            %center the window
            set( handles.figure2, 'Position', position );
            if handles.geometry ~= 0 ; copyobj(handles.geometry,handles.figure2);end
            
        end
        
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ ] = graphics_S23(io, P, t, Center, normals);
    %%   Program I/O graphics
    %    variable io is zero for the input graphics and one for the output graphics 
    xscale = max(P(:, 1))-min(P(:, 1));
    yscale = max(P(:, 2))-min(P(:, 2));
    zscale = max(P(:, 3))-min(P(:, 3));
    
    %%   Prepare common variables
    Ptemp = P'; ttemp = t'; ind = size(t, 1);  IND = 10000;
    X = reshape(Ptemp(1, ttemp(1:3, :)),[3, size(ttemp, 2)]);
    Y = reshape(Ptemp(2, ttemp(1:3, :)),[3, size(ttemp, 2)]);
    Z = reshape(Ptemp(3, ttemp(1:3, :)),[3, size(ttemp, 2)]);      
   
    %%   Input graphics (within the main GUI window)    
    if io>=0
        colormap jet;
        if ind<IND ec='y'; else ec='none'; end
        patch('Vertices', P, 'Faces', t(:, 1:3), 'FaceColor',  'b', 'FaceAlpha', 1.00, 'EdgeColor', ec)     
        xlabel('x, m'), ylabel('y, m'), zlabel('z, m'); 
        title('Problem geometry: object mesh');
        if zscale<max([xscale yscale])1e-6;
            temp = P;
            temp(:, 3) = 0.01*max([xscale yscale]);
            patch('Vertices', temp, 'Faces', t(1, 1:3), 'FaceColor',  'none', 'FaceAlpha', 0, 'EdgeColor', 'none')    
        end
            
    end
    if io>0
        hold on;
        quiver3(Center(:,1),Center(:,2),Center(:,3), ...
             normals(:,1),normals(:,2),normals(:,3), 0.5*sqrt(size(t, 1)/500), 'color', 'r');
    end

    
    %%   General settings 
    axis('equal'); axis('tight'); grid on; 
    view(-47, 20);
end

function [] = box(L, W, H, XC, YC, ZC, Color, Transparency, LineWidth, EdgeColor)
%%  Draws a rectangular object in 3D (a plane, a brick, or a line)
%   Usage: box(2, 2, 2, 0, 0, 0, [1 1 0], 0, 1, [ 0 0 0]);
%   SNM Summer 2012
%   __________L__________
%   |         |y        |
%   |W        |         |
%   |         *(XC,YC,ZC)   
%   |         |         |
%   |_________|_________|x
%       

    if (L>0) & (W>0)
        hr = patch([-L/2 -L/2 +L/2 +L/2]+XC, [-W/2 +W/2 +W/2 -W/2]+YC, [-H/2 -H/2 -H/2 -H/2]+ZC, ...
            Color, 'FaceAlpha', Transparency, 'LineWidth', LineWidth, 'EdgeColor', EdgeColor);   %   bottom
        hr = patch([-L/2 -L/2 +L/2 +L/2]+XC, [-W/2 +W/2 +W/2 -W/2]+YC, [+H/2 +H/2 +H/2 +H/2]+ZC, ...
            Color, 'FaceAlpha', Transparency, 'LineWidth', LineWidth, 'EdgeColor', EdgeColor);   %   top
    end
    if (W>0) & (H>0)    
        hr = patch([+L/2 +L/2 +L/2 +L/2]+XC, [-W/2 -W/2 +W/2 +W/2]+YC, [-H/2 +H/2 +H/2 -H/2]+ZC, ...
            Color, 'FaceAlpha', Transparency, 'LineWidth', LineWidth, 'EdgeColor', EdgeColor);   %   right   
        hr = patch([-L/2 -L/2 -L/2 -L/2]+XC, [-W/2 -W/2 +W/2 +W/2]+YC, [-H/2 +H/2 +H/2 -H/2]+ZC, ...
            Color, 'FaceAlpha', Transparency, 'LineWidth', LineWidth, 'EdgeColor', EdgeColor);   %   left
    end
    if (L>0) & (H>0)
        hr = patch([-L/2 -L/2 +L/2 +L/2]+XC, [+W/2 +W/2 +W/2 +W/2]+YC, [-H/2 +H/2 +H/2 -H/2]+ZC, ...
            Color, 'FaceAlpha', Transparency, 'LineWidth', LineWidth, 'EdgeColor', EdgeColor);   %   back   
        hr = patch([-L/2 -L/2 +L/2 +L/2]+XC, [-W/2 -W/2 -W/2 -W/2]+YC, [-H/2 +H/2 +H/2 -H/2]+ZC, ...
            Color, 'FaceAlpha', Transparency, 'LineWidth', LineWidth, 'EdgeColor', EdgeColor);   %   front
    end
    if (W==0)&(H==0)    %  x-line
        hrl = line([-L/2+XC L/2+XC], [YC YC], [ZC ZC], 'LineWidth', 3, 'Color', Color);  
    end
    if (L==0)&(H==0)    %  y-line
        hrl = line([XC XC], [-W/2+YC W/2+YC], [ZC ZC], 'LineWidth', 3, 'Color', Color);  
    end
    if (L==0)&(W==0)    %  z-line
        hrl = line([XC XC], [YC YC], [-H/2+ZC H/2+ZC], 'LineWidth', 3, 'Color', Color);  
    end
    h = 1; 

end

function [P, t] = brick(L, W, H, Tr, par)
    %   Create a uniform or non-uniform mesh (ratio 1:5) for a 
    %   base brick with the length L, width W, and height H 
    %   Tri is the approximate number of triangles for top/bottom plates
    %   SNM Summer 2012

    %   Generate the plate mesh first
    [P0, t0] = plate(L, W, Tr, par, 0, 0, 1);
    P0(:, 3) = [];
    
    %   Apply extrusion
    [P, t] = quadric(H, P0, t0); 
end

function [P, t] = circle(R, Tr, par, nx, ny, nz);
    %   Create a uniform or nonuniform (ratio 1:5) mesh for a base circle 
    %   of radius R, normal vector nx, ny, nz, and with approximately M 
    %   triangles
    %   Original code: DISTMESH 2004-2012 Per-Olof Persson
    %   Adopted and simplified: SNM Summer 2012 

    h    = waitbar(0.5, 'Please wait - computing triangular mesh'); 
    
    if Tr<20
        %   Create a simple structured triangular mesh
        M = Tr;
        x = R*[cos(2*pi*[0:M]/(M+1)) 0];
        y = R*[sin(2*pi*[0:M]/(M+1)) 0];
        P(:, 1) = x';
        P(:, 2) = y';
        P(:, 3) = 0;
        t(:, 1) = [1:M+1  ]';
        t(:, 2) = [2:M+1 1]';
        t(:, 3) = (M+2);
        close(h);
        return;
    end
    
    h0    = 2.5*R/sqrt(Tr);                             %   Approximate desired edge length
    xmin  = -R; xmax = +R; ymin = -R; ymax = +R;        %   Bounding box

    dptol   = 0.001;                                    %   Termination criterion
    ttol = 0.1;                                         %   Criterion for retriangulation                  
    Fscale  = 1.2;                                      %   Inner "internal pressure" (increased by 0.01)
    deltat  = 0.2;                                      %   Ratio between force and movement
    geps    = 0.001*h0;                                 %   Relative boundary tolerance
    deps    = sqrt(eps)*h0;                             %   Accuracy of gradient calculation

    %   Create initial distribution in bounding box 
    %   Equilateral triangles - faster convergence
    [x, y] = meshgrid(xmin:h0:xmax, ymin:h0*sqrt(3)/2:ymax);
    x(2:2:end, :) = x(2:2:end, :) + h0/2;               %   Shift even rows
    P = [x(:), y(:)];                                   %   List of vertices

    %   Remove vertices outside the region
    dist = dcircle(P, R);                               %   A column of signed distances
    P    = P(dist<geps, :);                             %   Keep only d<0 (geps) points
    N    = size(P, 1);                                  %   Number of vertices N
    
    Pold = inf;                                         %   For the first iteration
    count = 0;
    
    while 1        
        %   Retriangulation by the Delaunay algorithm
        if max(sqrt(sum((P-Pold).^2,2))/h0)>ttol        %   Any large movement?
            Pold    = P;                                %   Save current positions
            dt      = DelaunayTri(P);                   %   2D Delaunay triangulation
            t       = dt.Triangulation;                 %   List of triangles
            ic      = incenters(dt);                    %   Centroids of triangles
            dist    = dcircle(ic, R);                   %   A column of signed distances
            t       = t(dist<-geps, :);                 %   Keep interior trianges
            bars    = edges(TriRep(t, P));              %   All sorted mesh edges [M, 2]
            M       = size(bars, 1);                    %   Length of the array of edges       
        end
        count   = count + 1;

        %   Move mesh points based on edge lengths L1 and forces F
        barvec  = P(bars(:, 1),:)-P(bars(:, 2),:);            %   List of edge vectors [M, 2]
        L1       = sqrt(sum(barvec.^2,2));                    %   Edge lengths [M, 1]
        barsc   = 0.5*(P(bars(:, 1), :) + P(bars(:, 2), :));  %   Edge centers [M, 2]
        hbars   = hcircle(barsc, R, par);                     %   Element sizes for edge centers [M, 1]
        hbars   = hbars/sqrt(sum(hbars.^2)/M);                %   Element sizes normalized by its average [M, 1] 
        L0      = Fscale*sqrt(sum(L1.^2)/M)*hbars;            %   Desired (non-uniform) edge lengths [M, 1] 
        F       = max(L0-L1, 0);                              %   Edge forces [M, 1]
        Fvec    = repmat(F./L1, 1, 2).*barvec;                %   Edge forces (x,y components) [M, 2]

        Ftot = full(sparse(bars(:, [1,1,2,2]), ones(size(F))*[1,2,1,2], [Fvec, -Fvec], N, 2));
        P = P + deltat*Ftot;                                  %   Update node positions

        %   Bring outside vertices back to the boundary
        dist = dcircle(P, R);                           %   Find distances
        ind   = dist>0;                                 %   Find vertices outside (d>0)
        dgradx = (dcircle([P(ind, 1)+deps, P(ind, 2)], R)-dist(ind))/deps;   %   Numerical
        dgrady = (dcircle([P(ind, 1), P(ind, 2)+deps], R)-dist(ind))/deps;   %   gradient
        dgrad2 = dgradx.^2 + dgrady.^2;
        P(ind, :) = P(ind, :) - ...
            [dist(ind).*dgradx./dgrad2, dist(ind).*dgrady./dgrad2];     %   Projection

        %   Termination criterion: All interior nodes move less than dptol (scaled)
        factor = max(sqrt(sum(deltat*Ftot(dist<-geps,:).^2, 2))/h0); 
        if factor<dptol, break; end
        if count>2^12 break; end;
    end
    
    close(h);
    
       %   Rotate mesh
    P(:, 3) = 0;
    theta   = acos(nz)*sign(ny + eps);           %   rotation about the x-axis following the right-hand rule
    phi     = asin(nx/sqrt(nx^2+ny^2+eps^2));    %   rotation about the z-axis (the same way) 
    theta   = theta/pi*180;
    phi     = phi/pi*180; 
    P       = rotatex(P, theta);
    P       = rotatez(P, phi);
end
 
function [P, t] = combine(varargin)
    %%   Combine multiple meshes into one mesh 
    %   Usage: [P, t] = combine(P1, P2, P3, t1, t2, t3);
    %   SNM Summer 2012    
    MN = size(varargin, 2); %   total number of arguments
    
    %   Combine arrays of vertices
    Ptemp = [];
    for m = 1:MN/2
        Ptemp = [Ptemp varargin{m}'];
    end
    P = Ptemp';
   
    %   Combine arrays of triangles
    temp = 0;
    ttemp = [];
    for m = MN/2+1:MN        
        ttemp = [ttemp (varargin{m}'+temp)];
        temp = temp + size(varargin{m-MN/2}, 1);
    end
    t = ttemp';    
end

function [P, t] = cylinder(R, H, Tr, par);
    %   Create a uniform or non-uniform (ratio 1:5) mesh for a base 
    %   cylinder with the radius R, height h, and orientation about 
    %   the z-axis. Tri is approximate number of triangles
    %   SNM Summer 2012

     %   Generate the circle mesh first
    [P0, t0] = circle(R, Tr, par, 0, 0, 1);
    P0(:, 3) = [];
    
    %   Apply extrusion
    [P, t] = quadric(H, P0, t0); 
 
end

function d = dcircle(P, R)
    %   Compute signed distance function for a circle with radius R 
    %   Original code: DISTMESH 2004-2012 Per-Olof Persson (PhD thesis, pp. 38-39)

    d = sqrt(P(:, 1).^2 + P(:, 2).^2) - R;
end

function d = drectangle(P, L, W)
    %   Compute signed distance function for a rectangle with length L and
    %   width W
    %   Original code: DISTMESH 2004-2012 Per-Olof Persson

    x1 = -L/2;  x2 = +L/2;
    y1 = -W/2;  y2 = +W/2;
    d1 = y1 - P(:, 2);
    d2 =-y2 + P(:, 2);
    d3 = x1 - P(:, 1);
    d4 =-x2 + P(:, 1);

    d5 = sqrt(d1.^2 + d3.^2);
    d6 = sqrt(d1.^2 + d4.^2);
    d7 = sqrt(d2.^2 + d3.^2);
    d8 = sqrt(d2.^2 + d4.^2);

    d = -min(min(min(-d1, -d2), -d3), -d4);

    ix = d1>0 & d3>0;
    d(ix) = d5(ix);
    ix = d1>0 & d4>0;
    d(ix) = d6(ix);
    ix = d2>0 & d3>0;
    d(ix) = d7(ix);
    ix = d2>0 & d4>0;
    d(ix) = d8(ix);
end

function d = dsphere(P, R)
    %   Compute signed distance function for a sphere with radius R 
    %   Original code: DISTMESH 2004-2012 Per-Olof Persson (PhD thesis, pp. 38-39)

    d = sqrt(P(:, 1).^2 + P(:, 2).^2 + P(:, 3).^2) - R;
end

function h = hcircle(P, R, par)
    %   Compute element size function for a circle with radius R 
    %   0<par<1
    %   Original code: DISTMESH 2004-2012 Per-Olof Persson (PhD thesis, pp. 38-39)
    h = 1/par - (1/par-1)*sqrt(sum(P.^2, 2))/R;     %   Relative edge length
end

function h = hrectangle(P, L, W, par)
    %   Compute element size function for a rectangle
    %   0<par<1
    %   Original code: DISTMESH 2004-2012 Per-Olof Persson

    d = drectangle(P, L, W);
    h = 1/par - (1/par-1)*(L/2 + d).*(W/2 + d)/(L*W/4);    %   Relative edge length
end

function [Pnew] = move(P, x, y, z)
    %%   Move the center of the mesh by x y z
    %   SNM Summer 2012
    Pnew(:, 1) = +P(:, 1) + x;
    Pnew(:, 2) = +P(:, 2) + y;
    Pnew(:, 3) = +P(:, 3) + z;
end

function [centers, normals] = normcenters(P, t);
    %%   Find centers and outer normal vectors for convex bodies
    tr              = TriRep(t(:, 1:3), P);
    centers         = incenters(tr);
    normals         = faceNormals(tr);  %   normals from Tri class
    temp            = sign(dot(normals, centers, 2));
    temp(temp==0)   = 1;
    normals     = normals.*repmat(temp, 1, 3); 
end

function [P, t] = plate(L, W, Tr, par, nx, ny, nz)
    %   Create a uniform or nonuniform (ratio 1:5) mesh  for a base plate with 
    %   length L and width W, normal vector nx, ny, nz, and with approximately 
    %   Tr triangles
    %   Original code: DISTMESH 2004-2012 Per-Olof Persson
    %   Adopted and simplified: SNM Summer 2012 
    warning off;
    
    h    = waitbar(0.5, 'Please wait - computing triangular mesh');
    
    h0    = 1.40*sqrt(L*W)/sqrt(Tr);                    %   Approximate desired edge length
    xmin  = -L/2; xmax = +L/2; ymin = -W/2; ymax = +W/2;%   Bounding box

    dptol   = 0.001;                                    %   Termination criterion
    ttol    = 0.1;                                      %   Criterion for retriangulation  
    Fscale  = 1.2;                                      %   Inner "internal pressure" (original)
    if Tr<=200 & par<0.3
        Fscale = 1.10;
    end
    deltat  = 0.2;                                      %   Ratio between force and movement
    geps    = 0.001*h0;                                 %   Relative boundary tolerance (original)
    if par ==1
        geps    = 0.1*h0;                               %   Relative boundary tolerance (relaxed here)
    end  
    deps    = sqrt(eps)*h0;                             %   Accuracy of gradient calculation
    
    %   Create initial distribution in bounding box 
    %   Equilateral triangles - faster convergence
    [x, y] = meshgrid(xmin:h0:xmax, ymin:h0*sqrt(3)/2:ymax);
    x(2:2:end, :) = x(2:2:end, :) + h0/2;               %   Shift even rows
    P = [x(:), y(:)];                                   %   List of node coordinates

    %   Remove vertices outside the region
    dist = drectangle(P, L, W);                         %  A column of signed distances
    P    = P(dist<geps, :);                             %  Keep only d<0 (geps) points  
    N    = size(P, 1);                                  %  Number of vertices N
    Pold = inf;                                         %  For first iteration

    count   = 0;    
    while 1      
        %   Retriangulation by the Delaunay algorithm
        if max(sqrt(sum((P-Pold).^2,2))/h0)>ttol        %   Any large movement?
            Pold    = P;                                %   Save current positions
            dt      = DelaunayTri(P);                   %   2D Delaunay triangulation
            t       = dt.Triangulation;                 %   List of triangles
            ic      = incenters(dt);                    %   Centroids of triangles
            dist    = drectangle(ic, L, W);             %   A column of signed distances
            t       = t(dist<-geps, :);                 %   Keep interior trianges
            bars    = edges(TriRep(t, P));              %   All sorted mesh edges [M, 2]
            M       = size(bars, 1);                    %   Length of the array of edges        
        end
        count = count + 1;     

        %   Move mesh points based on edge lengths L and forces F
        barvec  = P(bars(:, 1),:)-P(bars(:, 2),:);              %   List of edge vectors [M, 2]
        L1      = sqrt(sum(barvec.^2,2));                       %   Edge lengths [M, 1]       
        barsc   = 0.5*(P(bars(:, 1), :) + P(bars(:, 2), :));    %   Edge centers [M, 2]
        hbars   = hrectangle(barsc, L, W, par);                 %   Element sizes for edge centers [M, 1]
        hbars   = hbars/sqrt(sum(hbars.^2)/M);                  %   Element sizes normalized by its average [M, 1] 
        L0      = Fscale*sqrt(sum(L1.^2)/M)*hbars;              %   Desired (non-uniform) edge lengths [M, 1]        
        F       = max(L0-L1, 0);                                %   Edge forces [M, 1]
        Fvec    = repmat(F./L1, 1, 2).*barvec;                  %   Edge forces (x,y components) [M, 2]

        Ftot = full(sparse(bars(:, [1,1,2,2]), ones(size(F))*[1,2,1,2], [Fvec, -Fvec], N, 2));    
        P = P + deltat*Ftot;                            %   Update node positions

        %   Bring outside vertices back to the boundary
        dist  = drectangle(P, L, W);                    %   Find distances
        ind   = dist>0;                                 %   Find vertices outside (d>0)
        dgradx = (drectangle([P(ind, 1)+deps, P(ind, 2)], L, W)-dist(ind))/deps;   %   Numerical
        dgrady = (drectangle([P(ind, 1), P(ind, 2)+deps], L, W)-dist(ind))/deps;   %   gradient
        dgrad2 = dgradx.^2 + dgrady.^2;
        P(ind, :) = P(ind, :) - ...
            [dist(ind).*dgradx./dgrad2, dist(ind).*dgrady./dgrad2];     %   Projection

        %   Termination criterion: All interior nodes move less than dptol (scaled)
        factor = max(sqrt(sum(deltat*Ftot(dist<-geps,:).^2, 2))/h0);
    
        if factor<dptol, break; end
        if count>2^13 break; end;
    end
    close(h);    
    
    %   Rotate mesh
    P(:, 3) = 0;
    theta   = acos(nz)*sign(ny + eps);           %   rotation about the x-axis following the right-hand rule
    phi     = asin(nx/sqrt(nx^2+ny^2+eps^2));    %   rotation about the z-axis (the same way) 
    theta   = theta/pi*180;
    phi     = phi/pi*180; 
    P       = rotatex(P, theta);
    P       = rotatez(P, phi);
end

function [P, t] = plate0(L, W, Nx, Ny, indicator, nx, ny, nz)
    %   Create structured mesh for a base plate on the size L by W with the normal vector nx, ny, nz  
    %   Use indicator=1 for a non-uniform mesh
    %   SNM Summer 2012
    %   Create mesh in the xy-plane
    if indicator == 0
        x = [0:Nx]/Nx*L;    %   uniform grid
        y = [0:Ny]/Ny*W;    %   uniform grid
    else
        x = L/2*(1 - cos(pi*[0:Nx]/Nx));    %   non-uniform grid
        y = W/2*(1 - cos(pi*[0:Ny]/Ny));    %   non-uniform grid
    end
    P(:, 1) = repmat(x, 1, length(y))';
    P(:, 2) = reshape(repmat(y, length(x), 1), 1, length(x)*length(y))';
    P(:, 1) = P(:, 1) - mean(P(:, 1));
    P(:, 2) = P(:, 2) - mean(P(:, 2));
    P(:, 3) = 0;
    %   Define triangles 
    t = []; temp = [];
    for n = 1:Ny
        for m = 1:Nx
            t1 = [m m+1    m+Nx+2]' + (n-1)*(Nx+1);
            t2 = [m m+Nx+1 m+Nx+2]' + (n-1)*(Nx+1);
            temp = [temp t1 t2];
        end
    end
    t = temp'; 
    %   Rotate mesh
    theta   = acos(nz)*sign(ny + eps);              %   rotation about the x-axis following the right-hand rule
    phi     = asin(nx/sqrt(nx^2+ny^2+eps));         %   rotation about the z-axis (the same way) 
    theta   = theta/pi*180;
    phi     = phi/pi*180; 
    P       = rotatex(P, theta);
    P       = rotatez(P, phi);
end

function [P, t] = quadric(H, P0, t0)
    H
    %   Create a cylinder (quadric) of height H from any planar mesh using
    %   extrusion
    %   SNM Summer 2012

    %   Find boundary vertices
    bars0     = freeBoundary(TriRep(t0, P0));
    barvec    = P0(bars0(:, 1), :)-P0(bars0(:, 2),: );    %   List of edge vectors [M, 2]
    Lavg      = mean(sqrt(sum(barvec.^2, 2)));            %   Average edge length
    Pb = P0(bars0', :); Pb(2:2:end, :) = [];              %   Boundary vertices
    pb = size(Pb, 1);   barsb = [1:pb; [2:pb 1]]';        %   Boundary connectivity
   
    %  Determine z-divisions   
    NN      = ceil(H/Lavg);                               %   Number of intervals along the height
    hsize   = H/NN;                                       %   Size of the z-divisions
    
    %   Replicate the cap mesh and its boundary; establish the connectivity
    P = P0; p0 = size(P0, 1); P(:, 3) = - H/2;  %   Bottom cap - vertices                                                                              
    t = t0;                                     %   Bottom cap - faces
    
    if NN==1              
        %   Connectivity: bottom to top
        t1(:, 1:2) = bars0;                     %   Lower nodes        
        t1(:, 3)   = bars0(:, 1) + p0;          %   Upper node
        t2(:, 1:2) = bars0       + p0;          %   Upper nodes
        t2(:, 3)   = bars0(:, 2);               %   Lower node
        t = [t' t1' t2']';
    else
        %   Connectivity: bottom to layer
        P(end+1:end+pb, 1:2) = Pb;
        P(end+1-pb:end, 3) = -H/2 + hsize;
        t1(:, 1:2) = bars0;                     %   Lower nodes        
        t1(:, 3)   = barsb(:, 1) + p0;          %   Upper node
        t2(:, 1:2) = barsb       + p0;          %   Upper nodes
        t2(:, 3)   = bars0(:, 2);               %   Lower node
        t = [t' t1' t2']';
        %   Connectivity: layer to layer
        for m =2:NN-1                                 
            P(end+1:end+pb, 1:2) = Pb;
            P(end+1-pb:end, 3) = -H/2 + hsize*m;
            t1(:, 1:2) = barsb       + p0 + pb*(m-2);   %   Lower nodes        
            t1(:, 3)   = barsb(:, 1) + p0 + pb*(m-1);   %   Upper node
            t2(:, 1:2) = barsb       + p0 + pb*(m-1);   %   Upper nodes
            t2(:, 3)   = barsb(:, 2) + p0 + pb*(m-2);   %   Lower node
            t = [t' t1' t2']';
        end
        %   Connectivity: layer to top
        t1(:, 1:2) = barsb       + p0 + pb*(NN-2);   %   Lower nodes        
        t1(:, 3)   = bars0(:, 1) + p0 + (NN-1)*pb;   %   Upper node
        t2(:, 1:2) = bars0       + p0 + (NN-1)*pb;   %   Upper nodes
        t2(:, 3)   = barsb(:, 2) + p0 + pb*(NN-2);   %   Lower node
        t = [t' t1' t2']';
    end
    %   Add the top cap
    P0(:, 3) = H/2; P = [P' P0']';
    t = [t' t0'+(NN-1)*pb+p0]';
end

function [Pnew] = rotatex(P, xangle)
    %%   Rotate the mesh about the x-axis
    %   The local coordinate system is used with the origin at the center of gravity
    %   SNM Summer 2012
    anglex = xangle/180*pi;
    %   Rotation about the x-axis  
    Pnew(:, 1) = +P(:, 1);
    Pnew(:, 2) = +P(:, 2)*cos(anglex) - P(:, 3)*sin(anglex);
    Pnew(:, 3) = +P(:, 2)*sin(anglex) + P(:, 3)*cos(anglex);
    Pnew(:, 2) = Pnew(:, 2);
    Pnew(:, 3) = Pnew(:, 3);
end

function [Pnew] = rotatey(P, yangle)
    %%   Rotate the mesh about the y-axis
    %   The local coordinate system is used with the origin at the center of gravity
    %   SNM Summer 2012
    angley = yangle/180*pi;
    %   Rotation about the y-axis    
    Pnew(:, 1) = +P(:, 1)*cos(angley) + P(:, 3)*sin(angley);
    Pnew(:, 2) = +P(:, 2);
    Pnew(:, 3) = -P(:, 1)*sin(angley) + P(:, 3)*cos(angley);
    Pnew(:, 1) = Pnew(:, 1);
    Pnew(:, 3) = Pnew(:, 3);
end

function [Pnew] = rotatez(P, zangle)
    %%   Rotate the mesh about the z-axis
    %   The local coordinate system is used with the origin at the center of gravity
    %   SNM Summer 2012
    anglez = zangle/180*pi;
    %   Rotation about the z-axis         
    Pnew(:, 1)  = +P(:, 1)*cos(anglez) - P(:, 2)*sin(anglez);
    Pnew(:, 2)  = +P(:, 1)*sin(anglez) + P(:, 2)*cos(anglez);
    Pnew(:, 3)  = +P(:, 3); 
    Pnew(:, 1) = Pnew(:, 1);
    Pnew(:, 2) = Pnew(:, 2);
end

function q = simpqual(P, t)
%   Facet quality - radius ratio 
%   Original code: DISTMESH 2004-2012 Per-Olof Persson

    a = sqrt(sum((P(t(:, 2), :) - P(t(:, 1), :)).^2,2));
    b = sqrt(sum((P(t(:, 3), :) - P(t(:, 1), :)).^2,2));
    c = sqrt(sum((P(t(:, 3), :) - P(t(:, 2), :)).^2,2));
    r = 1/2*sqrt((b+c-a).*(c+a-b).*(a+b-c)./(a+b+c));
    R = a.*b.*c./sqrt((a+b+c).*(b+c-a).*(c+a-b).*(a+b-c));
    q = 2*r./R;
    
end

function [P, t] = sphere(R, Tr);
    %   Create a uniform mesh for a base sphere of radius R with approximately M triangles
    %   Original code: DISTMESH 2004-2012 Per-Olof Persson
    %   Adopted and simplified: SNM Summer 2012 
   
    h    = waitbar(0.5, 'Please wait - computing triangular mesh');  
    
    h0    = 5.5*R/sqrt(Tr);                             %   Approximate desired edge length
    xmin  = -R; xmax = +R; ymin = -R; ymax = +R;        %   3D Bounding box
    zmin  = -R; zmax = +R; 

    dptol   = 0.005;                                    %   Termination criterion
    ttol    = 0.1;                                      %   Criterion for retriangulation
    Fscale  = 1.2;                                      %   Inner "internal pressure"
    deltat  = 0.2;                                      %   Ratio between force and movement
    geps    = 0.001*h0;                                 %   Relative boundary tolerance
    deps    = sqrt(eps)*h0;                             %   Accuracy of gradient calculation

    %   Create initial uniform distribution in bounding box (MATLAB isosurface from grid)
    [x, y, z] = ndgrid(xmin:h0:xmax, ymin:h0:ymax, zmin:h0:zmax);
    pv        = isosurface(x, y, z, reshape(dsphere([x(:), y(:), z(:)], R), size(x)), 0);
    P         = pv.vertices;
    N    = size(P, 1);                                  %  Number of vertices N
    Pold = inf;                                         %  For first iteration

    count = 0;
    while 1
        %   Retriangulation by the Delaunay algorithm
        if max(sqrt(sum((P-Pold).^2,2))/h0)>ttol        %   Any large movement?
            Pold    = P;                                %   Save current positions
            DT  = DelaunayTri(P);                       %   3D Delaunay triangulation
            t   = freeBoundary(DT);                     %   3D Triangular mesh
            bars    = edges(TriRep(t, P));              %   All sorted mesh edges [M, 2]
            M       = size(bars, 1);                    %   Length of the array of edges
            count = count + 1;  
        end

        %   Move mesh points based on edge lengths L1 and forces F
        barvec    = P(bars(:, 1),:)-P(bars(:, 2),:);            %   List of edge vectors [M, 2]
        L1         = sqrt(sum(barvec.^2,2));                    %   Edge lengths [M, 1]
        L0        = Fscale*sqrt(sum(L1.^2)/M)*ones(M, 1);       %   Desired (average) edge lengths [M, 1]
        F         = max(L0-L1, 0);                              %   Edge forces [M, 1]
        Fvec      = repmat(F./L1, 1, 3).*barvec;                %   Edge forces (x,y,z components) [M, 3]

        Ftot = full(sparse(bars(:,[1,1,1,2,2,2]), ones(size(F))*[1,2,3,1,2,3], [Fvec,-Fvec], N, 3));
        P = P + deltat*Ftot;                                    %   Update node positions

        %   Bring all vertices back to the boundary
        dist   = dsphere(P, R);                                 %   Find distances
        dgradx = (dsphere([P(:, 1)+deps, P(:, 2), P(:, 3)], R)-dist)/deps;   %   Numerical
        dgrady = (dsphere([P(:, 1), P(:, 2)+deps, P(:, 3)], R)-dist)/deps;   %   gradient
        dgradz = (dsphere([P(:, 1), P(:, 2), P(:, 3)+deps], R)-dist)/deps;   %   gradient
        dgrad2 = dgradx.^2 + dgrady.^2 + dgradz.^2;
        P = P - [dist.*dgradx./dgrad2, dist.*dgrady./dgrad2, dist.*dgradz./dgrad2];     
                                                                             %   Projection
        %   Termination criterion: All interior nodes move less than dptol (scaled)
        factor = max(sqrt(sum((P-Pold).^2, 2))/h0);
        if factor<dptol, break; end
        if count>2^10 break; end;
    end
    if exist('h') 
        close(h);
    end
end

