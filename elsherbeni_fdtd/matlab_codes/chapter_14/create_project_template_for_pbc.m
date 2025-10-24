function varargout = create_project_template_for_pbc(varargin)
% CREATE_PROJECT_TEMPLATE_FOR_PBC M-file for create_project_template_for_pbc.fig
%      CREATE_PROJECT_TEMPLATE_FOR_PBC, by itself, creates a new CREATE_PROJECT_TEMPLATE_FOR_PBC or raises the existing
%      singleton*.
%
%      H = CREATE_PROJECT_TEMPLATE_FOR_PBC returns the handle to a new CREATE_PROJECT_TEMPLATE_FOR_PBC or the handle to
%      the existing singleton*.
%
%      CREATE_PROJECT_TEMPLATE_FOR_PBC('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CREATE_PROJECT_TEMPLATE_FOR_PBC.M with the given input arguments.
%
%      CREATE_PROJECT_TEMPLATE_FOR_PBC('Property','Value',...) creates a new CREATE_PROJECT_TEMPLATE_FOR_PBC or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before create_project_template_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to create_project_template_for_pbc_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help create_project_template_for_pbc

% Last Modified by GUIDE v2.5 25-May-2012 11:08:48

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @create_project_template_for_pbc_OpeningFcn, ...
                   'gui_OutputFcn',  @create_project_template_for_pbc_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before create_project_template_for_pbc is made visible.
function create_project_template_for_pbc_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to create_project_template_for_pbc (see VARARGIN)

% Choose default command line output for create_project_template_for_pbc
handles.output = hObject;
%set(hObject,'color',[254 253 222]/255);
set(hObject,'color','w');
set(hObject,'name','Create project template for PBC');
addpath('fdtd_files');

set(handles.g_panel,'Parent',hObject);
set(handles.sl_panel,'Parent',hObject);
set(handles.o_panel,'Parent',hObject);

pos = get(hObject,'position');
pos(4) = 20;
set(hObject,'position',pos);

pos = get(handles.ps_panel,'position');
pos(2) = 3;
set(handles.ps_panel,'position',pos);

pos = get(handles.g_panel,'position');
pos(2) = 11;
set(handles.g_panel,'position',pos);

pos = get(handles.sl_panel,'position');
pos(2) = 8;
set(handles.sl_panel,'position',pos);

pos = get(handles.o_panel,'position');
pos(2) = 4.9;
set(handles.o_panel,'position',pos);

pos = get(handles.ps_toggle,'position');
pos(2) = 18;
set(handles.ps_toggle,'position',pos);

pos = get(handles.gsl_toggle,'position');
pos(2) = 18;
set(handles.gsl_toggle,'position',pos);

pos = get(handles.o_toggle,'position');
pos(2) = 18;
set(handles.o_toggle,'position',pos);

pos = get(handles.pname_text,'position');
pos(2) = 1;
set(handles.pname_text,'position',pos);

pos = get(handles.project_name,'position');
pos(2) = 1;
set(handles.project_name,'position',pos);

pos = get(handles.crbutton,'position');
pos(2) = 1;
set(handles.crbutton,'position',pos);

pos = get(handles.cbutton,'position');
pos(2) = 1;
set(handles.cbutton,'position',pos);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes create_project_template_for_pbc wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = create_project_template_for_pbc_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function number_of_time_steps_Callback(hObject, eventdata, handles)
% hObject    handle to number_of_time_steps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of number_of_time_steps as text
%        str2double(get(hObject,'String')) returns contents of number_of_time_steps as a double


% --- Executes during object creation, after setting all properties.
function number_of_time_steps_CreateFcn(hObject, eventdata, handles)
% hObject    handle to number_of_time_steps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function courant_factor_Callback(hObject, eventdata, handles)
% hObject    handle to courant_factor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of courant_factor as text
%        str2double(get(hObject,'String')) returns contents of courant_factor as a double


% --- Executes during object creation, after setting all properties.
function courant_factor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to courant_factor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function number_of_cells_per_wavelength_Callback(hObject, eventdata, handles)
% hObject    handle to number_of_cells_per_wavelength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of number_of_cells_per_wavelength as text
%        str2double(get(hObject,'String')) returns contents of number_of_cells_per_wavelength as a double


% --- Executes during object creation, after setting all properties.
function number_of_cells_per_wavelength_CreateFcn(hObject, eventdata, handles)
% hObject    handle to number_of_cells_per_wavelength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dx_Callback(hObject, eventdata, handles)
% hObject    handle to dx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dx as text
%        str2double(get(hObject,'String')) returns contents of dx as a double


% --- Executes during object creation, after setting all properties.
function dx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dy_Callback(hObject, eventdata, handles)
% hObject    handle to dy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dy as text
%        str2double(get(hObject,'String')) returns contents of dy as a double


% --- Executes during object creation, after setting all properties.
function dy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dz_Callback(hObject, eventdata, handles)
% hObject    handle to dz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dz as text
%        str2double(get(hObject,'String')) returns contents of dz as a double


% --- Executes during object creation, after setting all properties.
function dz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in type_xn.
function type_xn_Callback(hObject, eventdata, handles)
% hObject    handle to type_xn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns type_xn contents as cell array
%        contents{get(hObject,'Value')} returns selected item from type_xn


% --- Executes during object creation, after setting all properties.
function type_xn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to type_xn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function agnc_xn_Callback(hObject, eventdata, handles)
% hObject    handle to agnc_xn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of agnc_xn as text
%        str2double(get(hObject,'String')) returns contents of agnc_xn as a double


% --- Executes during object creation, after setting all properties.
function agnc_xn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to agnc_xn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cnc_xn_Callback(hObject, eventdata, handles)
% hObject    handle to cnc_xn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cnc_xn as text
%        str2double(get(hObject,'String')) returns contents of cnc_xn as a double


% --- Executes during object creation, after setting all properties.
function cnc_xn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cnc_xn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in type_yn.
function type_yn_Callback(hObject, eventdata, handles)
% hObject    handle to type_yn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns type_yn contents as cell array
%        contents{get(hObject,'Value')} returns selected item from type_yn


% --- Executes during object creation, after setting all properties.
function type_yn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to type_yn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function agnc_yn_Callback(hObject, eventdata, handles)
% hObject    handle to agnc_yn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of agnc_yn as text
%        str2double(get(hObject,'String')) returns contents of agnc_yn as a double


% --- Executes during object creation, after setting all properties.
function agnc_yn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to agnc_yn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cnc_yn_Callback(hObject, eventdata, handles)
% hObject    handle to cnc_yn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cnc_yn as text
%        str2double(get(hObject,'String')) returns contents of cnc_yn as a double


% --- Executes during object creation, after setting all properties.
function cnc_yn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cnc_yn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in type_zn.
function type_zn_Callback(hObject, eventdata, handles)
% hObject    handle to type_zn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns type_zn contents as cell array
%        contents{get(hObject,'Value')} returns selected item from type_zn


% --- Executes during object creation, after setting all properties.
function type_zn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to type_zn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function agnc_zn_Callback(hObject, eventdata, handles)
% hObject    handle to agnc_zn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of agnc_zn as text
%        str2double(get(hObject,'String')) returns contents of agnc_zn as a double


% --- Executes during object creation, after setting all properties.
function agnc_zn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to agnc_zn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cnc_zn_Callback(hObject, eventdata, handles)
% hObject    handle to cnc_zn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cnc_zn as text
%        str2double(get(hObject,'String')) returns contents of cnc_zn as a double


% --- Executes during object creation, after setting all properties.
function cnc_zn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cnc_zn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in type_xp.
function type_xp_Callback(hObject, eventdata, handles)
% hObject    handle to type_xp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns type_xp contents as cell array
%        contents{get(hObject,'Value')} returns selected item from type_xp


% --- Executes during object creation, after setting all properties.
function type_xp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to type_xp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function agnc_xp_Callback(hObject, eventdata, handles)
% hObject    handle to agnc_xp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of agnc_xp as text
%        str2double(get(hObject,'String')) returns contents of agnc_xp as a double


% --- Executes during object creation, after setting all properties.
function agnc_xp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to agnc_xp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cnc_xp_Callback(hObject, eventdata, handles)
% hObject    handle to cnc_xp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cnc_xp as text
%        str2double(get(hObject,'String')) returns contents of cnc_xp as a double


% --- Executes during object creation, after setting all properties.
function cnc_xp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cnc_xp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in type_yp.
function type_yp_Callback(hObject, eventdata, handles)
% hObject    handle to type_yp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns type_yp contents as cell array
%        contents{get(hObject,'Value')} returns selected item from type_yp


% --- Executes during object creation, after setting all properties.
function type_yp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to type_yp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function agnc_yp_Callback(hObject, eventdata, handles)
% hObject    handle to agnc_yp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of agnc_yp as text
%        str2double(get(hObject,'String')) returns contents of agnc_yp as a double


% --- Executes during object creation, after setting all properties.
function agnc_yp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to agnc_yp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cnc_yp_Callback(hObject, eventdata, handles)
% hObject    handle to cnc_yp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cnc_yp as text
%        str2double(get(hObject,'String')) returns contents of cnc_yp as a double


% --- Executes during object creation, after setting all properties.
function cnc_yp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cnc_yp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in type_zp.
function type_zp_Callback(hObject, eventdata, handles)
% hObject    handle to type_zp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns type_zp contents as cell array
%        contents{get(hObject,'Value')} returns selected item from type_zp


% --- Executes during object creation, after setting all properties.
function type_zp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to type_zp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function agnc_zp_Callback(hObject, eventdata, handles)
% hObject    handle to agnc_zp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of agnc_zp as text
%        str2double(get(hObject,'String')) returns contents of agnc_zp as a double


% --- Executes during object creation, after setting all properties.
function agnc_zp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to agnc_zp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cnc_zp_Callback(hObject, eventdata, handles)
% hObject    handle to cnc_zp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cnc_zp as text
%        str2double(get(hObject,'String')) returns contents of cnc_zp as a double


% --- Executes during object creation, after setting all properties.
function cnc_zp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cnc_zp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cpml_order_Callback(hObject, eventdata, handles)
% hObject    handle to cpml_order (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cpml_order as text
%        str2double(get(hObject,'String')) returns contents of cpml_order as a double


% --- Executes during object creation, after setting all properties.
function cpml_order_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cpml_order (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cpml_sigma_factor_Callback(hObject, eventdata, handles)
% hObject    handle to cpml_sigma_factor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cpml_sigma_factor as text
%        str2double(get(hObject,'String')) returns contents of cpml_sigma_factor as a double


% --- Executes during object creation, after setting all properties.
function cpml_sigma_factor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cpml_sigma_factor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cpml_kappa_max_Callback(hObject, eventdata, handles)
% hObject    handle to cpml_kappa_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cpml_kappa_max as text
%        str2double(get(hObject,'String')) returns contents of cpml_kappa_max as a double


% --- Executes during object creation, after setting all properties.
function cpml_kappa_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cpml_kappa_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cpml_alpha_min_Callback(hObject, eventdata, handles)
% hObject    handle to cpml_alpha_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cpml_alpha_min as text
%        str2double(get(hObject,'String')) returns contents of cpml_alpha_min as a double


% --- Executes during object creation, after setting all properties.
function cpml_alpha_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cpml_alpha_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cpml_alpha_max_Callback(hObject, eventdata, handles)
% hObject    handle to cpml_alpha_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cpml_alpha_max as text
%        str2double(get(hObject,'String')) returns contents of cpml_alpha_max as a double


% --- Executes during object creation, after setting all properties.
function cpml_alpha_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cpml_alpha_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in air.
function air_Callback(hObject, eventdata, handles)
% hObject    handle to air (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of air


% --- Executes on button press in pec.
function pec_Callback(hObject, eventdata, handles)
% hObject    handle to pec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of pec


% --- Executes on button press in pmc.
function pmc_Callback(hObject, eventdata, handles)
% hObject    handle to pmc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of pmc



function number_of_materials_Callback(hObject, eventdata, handles)
% hObject    handle to number_of_materials (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of number_of_materials as text
%        str2double(get(hObject,'String')) returns contents of number_of_materials as a double


% --- Executes during object creation, after setting all properties.
function number_of_materials_CreateFcn(hObject, eventdata, handles)
% hObject    handle to number_of_materials (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function number_of_bricks_Callback(hObject, eventdata, handles)
% hObject    handle to number_of_bricks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of number_of_bricks as text
%        str2double(get(hObject,'String')) returns contents of number_of_bricks as a double


% --- Executes during object creation, after setting all properties.
function number_of_bricks_CreateFcn(hObject, eventdata, handles)
% hObject    handle to number_of_bricks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function number_of_spheres_Callback(hObject, eventdata, handles)
% hObject    handle to number_of_spheres (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of number_of_spheres as text
%        str2double(get(hObject,'String')) returns contents of number_of_spheres as a double


% --- Executes during object creation, after setting all properties.
function number_of_spheres_CreateFcn(hObject, eventdata, handles)
% hObject    handle to number_of_spheres (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function number_of_thin_wires_Callback(hObject, eventdata, handles)
% hObject    handle to number_of_thin_wires (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of number_of_thin_wires as text
%        str2double(get(hObject,'String')) returns contents of number_of_thin_wires as a double


% --- Executes during object creation, after setting all properties.
function number_of_thin_wires_CreateFcn(hObject, eventdata, handles)
% hObject    handle to number_of_thin_wires (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function kx_start_Callback(hObject, eventdata, handles)
% hObject    handle to kx_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of kx_start as text
%        str2double(get(hObject,'String')) returns contents of kx_start as a double


% --- Executes during object creation, after setting all properties.
function kx_start_CreateFcn(hObject, eventdata, handles)
% hObject    handle to kx_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ky_start_Callback(hObject, eventdata, handles)
% hObject    handle to ky_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ky_start as text
%        str2double(get(hObject,'String')) returns contents of ky_start as a double


% --- Executes during object creation, after setting all properties.
function ky_start_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ky_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function number_of_resistors_Callback(hObject, eventdata, handles)
% hObject    handle to number_of_resistors (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of number_of_resistors as text
%        str2double(get(hObject,'String')) returns contents of number_of_resistors as a double


% --- Executes during object creation, after setting all properties.
function number_of_resistors_CreateFcn(hObject, eventdata, handles)
% hObject    handle to number_of_resistors (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function number_of_capacitors_Callback(hObject, eventdata, handles)
% hObject    handle to number_of_capacitors (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of number_of_capacitors as text
%        str2double(get(hObject,'String')) returns contents of number_of_capacitors as a double


% --- Executes during object creation, after setting all properties.
function number_of_capacitors_CreateFcn(hObject, eventdata, handles)
% hObject    handle to number_of_capacitors (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function number_of_inductors_Callback(hObject, eventdata, handles)
% hObject    handle to number_of_inductors (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of number_of_inductors as text
%        str2double(get(hObject,'String')) returns contents of number_of_inductors as a double


% --- Executes during object creation, after setting all properties.
function number_of_inductors_CreateFcn(hObject, eventdata, handles)
% hObject    handle to number_of_inductors (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function number_of_diodes_Callback(hObject, eventdata, handles)
% hObject    handle to number_of_diodes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of number_of_diodes as text
%        str2double(get(hObject,'String')) returns contents of number_of_diodes as a double


% --- Executes during object creation, after setting all properties.
function number_of_diodes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to number_of_diodes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in incident_plane_wave.
function incident_plane_wave_Callback(hObject, eventdata, handles)
% hObject    handle to incident_plane_wave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of incident_plane_wave



function number_of_electric_fields_Callback(hObject, eventdata, handles)
% hObject    handle to number_of_electric_fields (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of number_of_electric_fields as text
%        str2double(get(hObject,'String')) returns contents of number_of_electric_fields as a double


% --- Executes during object creation, after setting all properties.
function number_of_electric_fields_CreateFcn(hObject, eventdata, handles)
% hObject    handle to number_of_electric_fields (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function number_of_magnetic_fields_Callback(hObject, eventdata, handles)
% hObject    handle to number_of_magnetic_fields (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of number_of_magnetic_fields as text
%        str2double(get(hObject,'String')) returns contents of number_of_magnetic_fields as a double


% --- Executes during object creation, after setting all properties.
function number_of_magnetic_fields_CreateFcn(hObject, eventdata, handles)
% hObject    handle to number_of_magnetic_fields (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit38_Callback(hObject, eventdata, handles)
% hObject    handle to edit38 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit38 as text
%        str2double(get(hObject,'String')) returns contents of edit38 as a double


% --- Executes during object creation, after setting all properties.
function edit38_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit38 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit39_Callback(hObject, eventdata, handles)
% hObject    handle to edit39 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit39 as text
%        str2double(get(hObject,'String')) returns contents of edit39 as a double


% --- Executes during object creation, after setting all properties.
function edit39_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit39 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit40_Callback(hObject, eventdata, handles)
% hObject    handle to edit40 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit40 as text
%        str2double(get(hObject,'String')) returns contents of edit40 as a double


% --- Executes during object creation, after setting all properties.
function edit40_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit40 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function number_of_ports_Callback(hObject, eventdata, handles)
% hObject    handle to number_of_ports (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of number_of_ports as text
%        str2double(get(hObject,'String')) returns contents of number_of_ports as a double


% --- Executes during object creation, after setting all properties.
function number_of_ports_CreateFcn(hObject, eventdata, handles)
% hObject    handle to number_of_ports (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function number_of_sampled_voltages_Callback(hObject, eventdata, handles)
% hObject    handle to number_of_sampled_voltages (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of number_of_sampled_voltages as text
%        str2double(get(hObject,'String')) returns contents of number_of_sampled_voltages as a double


% --- Executes during object creation, after setting all properties.
function number_of_sampled_voltages_CreateFcn(hObject, eventdata, handles)
% hObject    handle to number_of_sampled_voltages (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function number_of_sampled_currents_Callback(hObject, eventdata, handles)
% hObject    handle to number_of_sampled_currents (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of number_of_sampled_currents as text
%        str2double(get(hObject,'String')) returns contents of number_of_sampled_currents as a double


% --- Executes during object creation, after setting all properties.
function number_of_sampled_currents_CreateFcn(hObject, eventdata, handles)
% hObject    handle to number_of_sampled_currents (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function frequency_start_Callback(hObject, eventdata, handles)
% hObject    handle to frequency_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of frequency_start as text
%        str2double(get(hObject,'String')) returns contents of frequency_start as a double


% --- Executes during object creation, after setting all properties.
function frequency_start_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frequency_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function frequency_step_Callback(hObject, eventdata, handles)
% hObject    handle to frequency_step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of frequency_step as text
%        str2double(get(hObject,'String')) returns contents of frequency_step as a double


% --- Executes during object creation, after setting all properties.
function frequency_step_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frequency_step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function frequency_end_Callback(hObject, eventdata, handles)
% hObject    handle to frequency_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of frequency_end as text
%        str2double(get(hObject,'String')) returns contents of frequency_end as a double


% --- Executes during object creation, after setting all properties.
function frequency_end_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frequency_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in run_simulation.
function run_simulation_Callback(hObject, eventdata, handles)
% hObject    handle to run_simulation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of run_simulation


% --- Executes on button press in show_material_mesh.
function show_material_mesh_Callback(hObject, eventdata, handles)
% hObject    handle to show_material_mesh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of show_material_mesh


% --- Executes on button press in show_problem_space.
function show_problem_space_Callback(hObject, eventdata, handles)
% hObject    handle to show_problem_space (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of show_problem_space



function plotting_step_Callback(hObject, eventdata, handles)
% hObject    handle to plotting_step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of plotting_step as text
%        str2double(get(hObject,'String')) returns contents of plotting_step as a double


% --- Executes during object creation, after setting all properties.
function plotting_step_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plotting_step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ps_labels.
function ps_labels_Callback(hObject, eventdata, handles)
% hObject    handle to ps_labels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ps_labels


% --- Executes on button press in ps_axis_at_origin.
function ps_axis_at_origin_Callback(hObject, eventdata, handles)
% hObject    handle to ps_axis_at_origin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ps_axis_at_origin


% --- Executes on button press in ps_axis_outside.
function ps_axis_outside_Callback(hObject, eventdata, handles)
% hObject    handle to ps_axis_outside (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ps_axis_outside


% --- Executes on button press in ps_outer_boundaries.
function ps_outer_boundaries_Callback(hObject, eventdata, handles)
% hObject    handle to ps_outer_boundaries (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ps_outer_boundaries


% --- Executes on button press in ps_cpml_boundaries.
function ps_cpml_boundaries_Callback(hObject, eventdata, handles)
% hObject    handle to ps_cpml_boundaries (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ps_cpml_boundaries


% --- Executes on button press in ps_xn_grid.
function ps_xn_grid_Callback(hObject, eventdata, handles)
% hObject    handle to ps_xn_grid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ps_xn_grid


% --- Executes on button press in ps_yn_grid.
function ps_yn_grid_Callback(hObject, eventdata, handles)
% hObject    handle to ps_yn_grid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ps_yn_grid


% --- Executes on button press in ps_zn_grid.
function ps_zn_grid_Callback(hObject, eventdata, handles)
% hObject    handle to ps_zn_grid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ps_zn_grid


% --- Executes on button press in ps_xp_grid.
function ps_xp_grid_Callback(hObject, eventdata, handles)
% hObject    handle to ps_xp_grid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ps_xp_grid


% --- Executes on button press in ps_yp_grid.
function ps_yp_grid_Callback(hObject, eventdata, handles)
% hObject    handle to ps_yp_grid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ps_yp_grid


% --- Executes on button press in ps_zp_grid.
function ps_zp_grid_Callback(hObject, eventdata, handles)
% hObject    handle to ps_zp_grid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ps_zp_grid



function number_of_animations_Callback(hObject, eventdata, handles)
% hObject    handle to number_of_animations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of number_of_animations as text
%        str2double(get(hObject,'String')) returns contents of number_of_animations as a double


% --- Executes during object creation, after setting all properties.
function number_of_animations_CreateFcn(hObject, eventdata, handles)
% hObject    handle to number_of_animations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)

nsrc = str2num(get(handles.kx_start,'string'));

if ~isfield(handles, 'vs_waveforms')
    handles.vs_waveforms = 3*ones(1, nsrc);
else
    sz = size(handles.vs_waveforms, 2);
    if sz < nsrc
        for mi=sz+1: number_of_sources
            handles.vs_waveforms(mi) = 3;
        end
    end
    if sz > nsrc
        handles.vs_waveforms(nsrc+1:sz) = [];
    end
end
ret = source_waveforms(nsrc, handles.vs_waveforms);
if ret~=0
    handles.vs_waveforms = ret;
end
guidata(hObject, handles);

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
nsrc = str2num(get(handles.ky_start,'string'));

if ~isfield(handles, 'cs_waveforms')
    handles.cs_waveforms = 3*ones(1, nsrc);
else
    sz = size(handles.cs_waveforms, 2);
    if sz < nsrc
        for mi=sz+1: nsrc
            handles.cs_waveforms(mi) = 3;
        end
    end
    if sz > nsrc
        handles.cs_waveforms(nsrc+1:sz) = [];
    end
end
ret = source_waveforms(nsrc, handles.cs_waveforms);
if ret~=0
    handles.cs_waveforms = ret;
end
guidata(hObject, handles);

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)

if ~isfield(handles, 'in_waveform')
    handles.in_waveform = 3;
end
ret = source_waveforms(1, handles.in_waveform);
if ret~=0
    handles.in_waveform = ret;
end
guidata(hObject, handles);

% --- Executes on button press in checkbox19.
function checkbox19_Callback(hObject, eventdata, handles)



% --- Executes on button press in cbutton.
function cbutton_Callback(hObject, eventdata, handles)

eval(get(gcf,'CloseRequestFcn'));

% --- Executes on button press in crbutton.
function crbutton_Callback(hObject, eventdata, handles)

project_name = get(handles.project_name,'string');
if exist(project_name,'dir')
    button = questdlg('The project directory exists. Do you want to overwrite?');
    if ~strcmp(button,'Yes')
        return;
    end
end
ret= mkdir(project_name);
if ~ret
    msgbox('An error occured. Can not create the project directory!');
    return;
end
cd(project_name);
handles = create_problem_space_file(handles);
handles = create_geometry_file(handles);
handles = create_sources_file(handles);
handles = create_outputs_file(handles);
cd('..');
% ================================================================
function handles = create_sources_file(handles)

fid = fopen('define_sources_and_lumped_elements.m','w');
if fid <= 0
    msgbox('An error occured. Can not create files!');
end

str = 'disp(''defining sources''); \n\n';
fprintf(fid,str);
  
ind = 0;


fclose(fid);

% ================================================================
% ================================================================
function handles = create_geometry_file(handles)

fid = fopen('define_geometry.m','w');
if fid <= 0
    msgbox('An error occured. Can not create files!');
end

str = 'disp(''defining the problem geometry''); \n\n';
fprintf(fid,str);

ind = 0;
nobj = str2num(get(handles.number_of_bricks,'string'));
if nobj>0
    for mi=1:nobj
        ind = ind+1; sind = num2str(ind);
        fprintf(fid,['%% a brick\n']);
        fprintf(fid,['bricks(' sind ').min_x = ' num2str(ind/10) ';\n']);
        fprintf(fid,['bricks(' sind ').min_y = ' num2str(ind/10) ';\n']);
        fprintf(fid,['bricks(' sind ').min_z = ' num2str(0) ';\n']);
        fprintf(fid,['bricks(' sind ').max_x = ' num2str(ind/10+0.1) ';\n']);
        fprintf(fid,['bricks(' sind ').max_y = ' num2str(ind/10+0.1) ';\n']);
        fprintf(fid,['bricks(' sind ').max_z = ' num2str(1/10) ';\n']);
        fprintf(fid,['bricks(' sind ').material_type = ' num2str(ind) ';\n\n']);
    end
end

ind = 0;
nobj = str2num(get(handles.number_of_spheres,'string'));
if nobj>0
    for mi=1:nobj
        ind = ind+1; sind = num2str(ind);
        fprintf(fid,['%% a sphere\n']);
        fprintf(fid,['spheres(' sind ').radius = ' num2str(ind/10) ';\n']);
        fprintf(fid,['spheres(' sind ').center_x = ' num2str(-1/10) ';\n']);
        fprintf(fid,['spheres(' sind ').center_y = ' num2str(-1/10) ';\n']);
        fprintf(fid,['spheres(' sind ').center_z = ' num2str(ind/5) ';\n']);
        fprintf(fid,['spheres(' sind ').material_type = ' num2str(ind) ';\n\n']);
    end
end

ind = 0;
nobj = str2num(get(handles.number_of_thin_wires,'string'));
if nobj>0
    for mi=1:nobj
        ind = ind+1; sind = num2str(ind);
        fprintf(fid,['%% a thin wire\n']);
        fprintf(fid,['thin_wires(' sind ').min_x = ' num2str(ind/10) ';\n']);
        fprintf(fid,['thin_wires(' sind ').min_y = ' num2str(-ind/10) ';\n']);
        fprintf(fid,['thin_wires(' sind ').min_z = ' num2str(0) ';\n']);
        fprintf(fid,['thin_wires(' sind ').max_x = ' num2str(ind/10) ';\n']);
        fprintf(fid,['thin_wires(' sind ').max_y = ' num2str(-ind/10) ';\n']);
        fprintf(fid,['thin_wires(' sind ').max_z = ' num2str(1/10) ';\n']);
        dx = str2num(get(handles.dx,'string'));
        fprintf(fid,['thin_wires(' sind ').radius = ' num2str(dx/4) ';\n']);
        fprintf(fid,['thin_wires(' sind ').direction = ''z'';\n\n']);
    end
end

fclose(fid);

% ================================================================
function handles = create_problem_space_file(handles)

fid = fopen('define_problem_space_parameters.m','w');
if fid <= 0
    msgbox('An error occured. Can not create files!');
end

str = 'disp(''defining the problem space parameters''); \n\n';
fprintf(fid,str);
str = '%% maximum number of time steps to run FDTD simulation \n';
fprintf(fid,str);
str = ['number_of_time_steps = ' get(handles.number_of_time_steps, 'string') ';\n\n']; 
fprintf(fid,str);

str = '%% A factor that determines duration of a time step wrt CFL limit\n';
fprintf(fid,str);
str = ['courant_factor = ' get(handles.courant_factor, 'string') ';\n\n']; 
fprintf(fid,str);

str = '%% Dimensions of a unit cell in x, y, and z directions (meters)\n';
fprintf(fid,str);
str = ['dx = ' get(handles.dx, 'string') ';\n'];    
fprintf(fid,str);
str = ['dy = ' get(handles.dy, 'string') ';\n'];    
fprintf(fid,str);
str = ['dz = ' get(handles.dz, 'string') ';\n\n'];    
fprintf(fid,str);


    str = '%% ==<periodic boundary simulation parameters>========\n';
    fprintf(fid,str);
    str = 'disp(''Defining periodic boundary simulation''); \n\n';
    fprintf(fid,str);
    
    if get(handles.TE_mode,'value') == 1
        fprintf(fid,'periodic_boundary.mode = ''TE'';\n');
        fprintf(fid,['periodic_boundary.E_phi = ' get(handles.E_phi,'string') ';\n']);
        fprintf(fid,['periodic_boundary.kx = ' get(handles.kx_start,'string') ';\n']);
        fprintf(fid,['periodic_boundary.ky = ' get(handles.ky_start,'string') ';\n']);
    end
    if get(handles.TM_mode,'value') == 1
        fprintf(fid,'periodic_boundary.mode = ''TM'';\n');
        fprintf(fid,['periodic_boundary.H_phi = ' get(handles.H_phi,'string') ';\n']);
        fprintf(fid,['periodic_boundary.kx = ' get(handles.kx_start,'string') ';\n']);
        fprintf(fid,['periodic_boundary.ky = ' get(handles.ky_start,'string') ';\n']);
    end
    if get(handles.TEM_mode,'value') == 1
        fprintf(fid,'periodic_boundary.mode = ''TEM'';\n');
        fprintf(fid,['periodic_boundary.E_x = ' get(handles.E_x,'string') ';\n']);
        fprintf(fid,['periodic_boundary.E_y = ' get(handles.E_y,'string') ';\n']);
    end
    fprintf(fid,['periodic_boundary.source_z = ' get(handles.source_z,'string') ';\n']);

    if get(handles.reflection_on,'value') == 1
        fprintf(fid,['periodic_boundary.reflection_z = ' get(handles.reflection_z,'string') ';\n']);
    end
        
    if get(handles.transmission_on,'value') == 1
        fprintf(fid,['periodic_boundary.transmission_z = ' get(handles.transmission_z,'string') ';\n']);
    end
    fprintf(fid,'\n');
        
str = '%% ==<boundary conditions>========\n';fprintf(fid,str);
str = '%% Here we define the boundary conditions parameters \n';fprintf(fid,str);
str = '%% ''pec'' : perfect electric conductor\n';fprintf(fid,str);
str = '%% ''cpml'' : convolutional PML \n';fprintf(fid,str);
str = '%% ''pbc'' : Periodic Boundary Condition \n';fprintf(fid,str);
str = '%% if cpml_number_of_cells is less than zero\n';fprintf(fid,str);
str = '%% CPML extends inside of the domain rather than outwards\n\n';fprintf(fid,str);

str = strvcat('cpml','pec','pbc');
st = strtrim(str(get(handles.type_xn,'value'),:));
fprintf(fid,['boundary.type_xn = ''' st ''';\n']);
fprintf(fid,['boundary.air_buffer_number_of_cells_xn = ' get(handles.agnc_xn,'string') ';\n']);
fprintf(fid,['boundary.cpml_number_of_cells_xn = ' get(handles.cnc_xn,'string') ';\n\n']);

st = strtrim(str(get(handles.type_xp,'value'),:));
fprintf(fid,['boundary.type_xp = ''' st ''';\n']);
fprintf(fid,['boundary.air_buffer_number_of_cells_xp = ' get(handles.agnc_xp,'string') ';\n']);
fprintf(fid,['boundary.cpml_number_of_cells_xp = ' get(handles.cnc_xp,'string') ';\n\n']);

st = strtrim(str(get(handles.type_yn,'value'),:));
fprintf(fid,['boundary.type_yn = ''' st ''';\n']);
fprintf(fid,['boundary.air_buffer_number_of_cells_yn = ' get(handles.agnc_yn,'string') ';\n']);
fprintf(fid,['boundary.cpml_number_of_cells_yn = ' get(handles.cnc_yn,'string') ';\n\n']);

st = strtrim(str(get(handles.type_yp,'value'),:));
fprintf(fid,['boundary.type_yp = ''' st ''';\n']);
fprintf(fid,['boundary.air_buffer_number_of_cells_yp = ' get(handles.agnc_yp,'string') ';\n']);
fprintf(fid,['boundary.cpml_number_of_cells_yp = ' get(handles.cnc_yp,'string') ';\n\n']);

st = strtrim(str(get(handles.type_zn,'value'),:));
fprintf(fid,['boundary.type_zn = ''' st ''';\n']);
fprintf(fid,['boundary.air_buffer_number_of_cells_zn = ' get(handles.agnc_zn,'string') ';\n']);
fprintf(fid,['boundary.cpml_number_of_cells_zn = ' get(handles.cnc_zn,'string') ';\n\n']);

st = strtrim(str(get(handles.type_zp,'value'),:));
fprintf(fid,['boundary.type_zp = ''' st ''';\n']);
fprintf(fid,['boundary.air_buffer_number_of_cells_zp = ' get(handles.agnc_zp,'string') ';\n']);
fprintf(fid,['boundary.cpml_number_of_cells_zp = ' get(handles.cnc_zp,'string') ';\n\n']);

fprintf(fid,['boundary.cpml_order = ' get(handles.cpml_order,'string') ';\n']);
fprintf(fid,['boundary.cpml_sigma_factor = ' get(handles.cpml_sigma_factor,'string') ';\n']);
fprintf(fid,['boundary.cpml_kappa_max = ' get(handles.cpml_kappa_max,'string') ';\n']);
fprintf(fid,['boundary.cpml_alpha_min = ' get(handles.cpml_alpha_min,'string') ';\n']);
fprintf(fid,['boundary.cpml_alpha_max = ' get(handles.cpml_alpha_max,'string') ';\n\n']);

fprintf(fid,'%% ===<material types>============\n');
fprintf(fid,'%% Here we define and initialize the arrays of material types\n');
fprintf(fid,'%% eps_r   : relative permittivity\n');
fprintf(fid,'%% mu_r    : relative permeability\n');
fprintf(fid,'%% sigma_e : electric conductivity\n');
fprintf(fid,'%% sigma_m : magnetic conductivity\n\n');

ind = 0;
if get(handles.air,'value')
    ind = ind+1; sind = num2str(ind);
    fprintf(fid,'%% air\n');
    fprintf(fid,['material_types(' sind ').eps_r   = 1;\n']);
    fprintf(fid,['material_types(' sind ').mu_r    = 1;\n']);
    fprintf(fid,['material_types(' sind ').sigma_e = 0;\n']);
    fprintf(fid,['material_types(' sind ').sigma_m = 0;\n']); 
    fprintf(fid,['material_types(' sind ').color   = [1 1 1];\n\n']);
end

if get(handles.pec,'value')
    ind = ind+1; sind = num2str(ind);
    fprintf(fid,'%% PEC : perfect electric conductor;\n');
    fprintf(fid,['material_types(' sind ').eps_r   = 1;\n']);
    fprintf(fid,['material_types(' sind ').mu_r    = 1;\n']);
    fprintf(fid,['material_types(' sind ').sigma_e = 1e10;\n']);
    fprintf(fid,['material_types(' sind ').sigma_m = 0;\n']); 
    fprintf(fid,['material_types(' sind ').color   = [1 0 0];\n\n']);
end

if get(handles.pmc,'value')
    ind = ind+1; sind = num2str(ind);
    fprintf(fid,['%% PMC : perfect magnetic conductor\n']);
    fprintf(fid,['material_types(' sind ').eps_r   = 1;\n']);
    fprintf(fid,['material_types(' sind ').mu_r    = 1;\n']);
    fprintf(fid,['material_types(' sind ').sigma_e = 0;\n']);
    fprintf(fid,['material_types(' sind ').sigma_m = 1e10;\n']);
    fprintf(fid,['material_types(' sind ').color   = [0 1 0];\n\n']);
end
nmat = str2num(get(handles.number_of_materials,'string'));
for mi=1:nmat
    ind = ind+1; sind = num2str(ind);
    fprintf(fid,['%% a material type\n']);
    fprintf(fid,['material_types(' sind ').eps_r   = ' num2str(ind) ';\n']);
    fprintf(fid,['material_types(' sind ').mu_r    = 1;\n']);
    fprintf(fid,['material_types(' sind ').sigma_e = 0;\n']);
    fprintf(fid,['material_types(' sind ').sigma_m = 1e10;\n']);
    r = num2str(mod(ind*80,255)/255);
    g = num2str(mod(ind*120,255)/255);
    b = num2str(mod(ind*160,255)/255);
    fprintf(fid,['material_types(' sind ').color   = [' r ' ' g ' ' b '];\n\n']);
end
fclose(fid);

% ================================================================
function handles = create_outputs_file(handles)

fid = fopen('define_output_parameters.m','w');
if fid <= 0
    msgbox('An error occured. Can not create files!');
end

str = 'disp(''defining output parameters''); \n\n';
fprintf(fid,str);


fprintf(fid,'%% figure refresh rate\n');
fprintf(fid,['plotting_step = ' get(handles.plotting_step,'string') ';\n\n']);

fprintf(fid,'%% mode of operation\n');
if get(handles.run_simulation,'value')
    fprintf(fid,'run_simulation = true;\n');
else
    fprintf(fid,'run_simulation = false;\n');
end
if get(handles.show_material_mesh,'value')
    fprintf(fid,'show_material_mesh = true;\n');
else
    fprintf(fid,'show_material_mesh = false;\n');
end
if get(handles.show_problem_space ,'value')
    fprintf(fid,'show_problem_space = true;\n');
else
    fprintf(fid,'show_problem_space = false;\n');
end
fprintf(fid,'\n');   

fprintf(fid,'%% frequency domain parameters\n');
fprintf(fid,['frequency_domain.start = ' get(handles.frequency_start,'string') ';\n']);
fprintf(fid,['frequency_domain.end   = ' get(handles.frequency_end,'string') ';\n']);
fprintf(fid,['frequency_domain.step  = ' get(handles.frequency_step,'string') ';\n']);


ind = 0;
nsrc = str2num(get(handles.number_of_electric_fields,'string'));
if nsrc>0
    fprintf(fid,'%% define sampled electric fields\n');
    fprintf(fid,'%% component: vector component ''x'',''y'',''z'', or magnitude ''m''\n');
    fprintf(fid,'%% display plot = true, in order to plot field during simulation \n\n');    
    
    for mi=1:nsrc
        ind = ind+1; sind = num2str(ind);
        fprintf(fid,['%% a sampled electric field\n']);
        fprintf(fid,['sampled_electric_fields(' sind ').x = ' num2str(ind/10) ';\n']);
        fprintf(fid,['sampled_electric_fields(' sind ').y = ' num2str(ind/10) ';\n']);
        fprintf(fid,['sampled_electric_fields(' sind ').z = ' num2str(0) ';\n']);
        fprintf(fid,['sampled_electric_fields(' sind ').component = ''z'';\n']);
        fprintf(fid,['sampled_electric_fields(' sind ').display_plot = true;\n']);
        fprintf(fid,'\n');
    end
end
fprintf(fid,'\n');
 
ind = 0;
nsrc = str2num(get(handles.number_of_magnetic_fields,'string'));
if nsrc>0
    fprintf(fid,'%% define sampled magnetic fields\n');
    fprintf(fid,'%% component: vector component ''x'',''y'',''z'', or magnitude ''m''\n');
    fprintf(fid,'%% display plot = true, in order to plot field during simulation \n\n');    
    
    for mi=1:nsrc
        ind = ind+1; sind = num2str(ind);
        fprintf(fid,['%% a sampled magnetic field\n']);
        fprintf(fid,['sampled_magnetic_fields(' sind ').x = ' num2str(ind/10) ';\n']);
        fprintf(fid,['sampled_magnetic_fields(' sind ').y = ' num2str(ind/10) ';\n']);
        fprintf(fid,['sampled_magnetic_fields(' sind ').z = ' num2str(0) ';\n']);
        fprintf(fid,['sampled_magnetic_fields(' sind ').component = ''z'';\n']);
        fprintf(fid,['sampled_magnetic_fields(' sind ').display_plot = true;\n']);
        fprintf(fid,'\n');
    end
end
fprintf(fid,'\n');


ind = 0;
nsrc = str2num(get(handles.number_of_sampled_voltages,'string'));
if nsrc>0    
    for mi=1:nsrc
        ind = ind+1; sind = num2str(ind);
        fprintf(fid,['%% a sampled voltage\n']);
        fprintf(fid,['sampled_voltages(' sind ').min_x = ' num2str(ind/10) ';\n']);
        fprintf(fid,['sampled_voltages(' sind ').min_y = ' num2str(ind/10) ';\n']);
        fprintf(fid,['sampled_voltages(' sind ').min_z = ' num2str(0) ';\n']);
        fprintf(fid,['sampled_voltages(' sind ').max_x = ' num2str(ind/10) ';\n']);
        fprintf(fid,['sampled_voltages(' sind ').max_y = ' num2str(ind/10) ';\n']);
        fprintf(fid,['sampled_voltages(' sind ').max_z = ' num2str(ind/10) ';\n']);
        fprintf(fid,['sampled_voltages(' sind ').direction = ''zp'';\n']);
        fprintf(fid,['sampled_voltages(' sind ').display_plot = true;\n']);
        fprintf(fid,'\n');
    end
end
fprintf(fid,'\n');

ind = 0;
nsrc = str2num(get(handles.number_of_sampled_currents,'string'));
if nsrc>0    
    for mi=1:nsrc
        ind = ind+1; sind = num2str(ind);
        fprintf(fid,['%% a sampled current\n']);
        fprintf(fid,['sampled_currents(' sind ').min_x = ' num2str(ind/10) ';\n']);
        fprintf(fid,['sampled_currents(' sind ').min_y = ' num2str(ind/10) ';\n']);
        fprintf(fid,['sampled_currents(' sind ').min_z = ' num2str(0) ';\n']);
        fprintf(fid,['sampled_currents(' sind ').max_x = ' num2str(ind/10) ';\n']);
        fprintf(fid,['sampled_currents(' sind ').max_y = ' num2str(ind/10) ';\n']);
        fprintf(fid,['sampled_currents(' sind ').max_z = ' num2str(ind/10) ';\n']);
        fprintf(fid,['sampled_currents(' sind ').direction = ''zp'';\n']);
        fprintf(fid,['sampled_currents(' sind ').display_plot = true;\n']);
        fprintf(fid,'\n');
    end
end
fprintf(fid,'\n');

ind = 0;
nsrc = str2num(get(handles.number_of_ports,'string'));
if nsrc>0    
    for mi=1:nsrc
        ind = ind+1; sind = num2str(ind);
        fprintf(fid,['%% a port\n']);
        fprintf(fid,['ports(' sind ').sampled_voltage_index = ' num2str(ind) ';\n']);
        fprintf(fid,['ports(' sind ').sampled_current_index = ' num2str(ind) ';\n']);
        fprintf(fid,['ports(' sind ').impedance = 50;\n']);
        if mi==1
            fprintf(fid,['ports(' sind ').is_source_port = true;\n']);
        else
            fprintf(fid,['ports(' sind ').is_source_port = false;\n']);
        end
        fprintf(fid,'\n');
    end
end
fprintf(fid,'\n');


ind = 0;
nsrc = str2num(get(handles.number_of_animations,'string'));
if nsrc>0    
fprintf(fid,'%% define animation\n');
fprintf(fid,'%% field_type shall be ''e'' or ''h''\n');
fprintf(fid,'%% plane cut shall be ''xy'', ''yz'', or ''zx''\n');
fprintf(fid,'%% component shall be ''x'', ''y'', ''z'', or ''m''\n');

    for mi=1:nsrc
        ind = ind+1; sind = num2str(ind);
        fprintf(fid,['%% an animation\n']);
        fprintf(fid,['animation(' sind ').field_type = ''e'';\n']);
        fprintf(fid,['animation(' sind ').component = ''m'';\n']);
        fprintf(fid,['animation(' sind ').plane_cut(1).type = ''xy'';\n']);
        fprintf(fid,['animation(' sind ').plane_cut(1).position  = 0;\n']);
        fprintf(fid,['animation(' sind ').display_grid = true;\n']);
        fprintf(fid,['animation(' sind ').display_objects = true;\n']);
        fprintf(fid,['animation(' sind ').save_movie = true;\n']);
        fprintf(fid,['animation(' sind ').view_angles = [40 30];\n']);
        fprintf(fid,['animation(' sind ').zoom = 0.8;\n']);
        fprintf(fid,['animation(' sind ').enable = true;\n']);
        fprintf(fid,'\n');
    end
end
fprintf(fid,'\n');

fprintf(fid,'%% display problem space parameters\n');
if get(handles.ps_labels,'value') 
    fprintf(fid,'problem_space_display.labels = true;\n');
else
    fprintf(fid,'problem_space_display.labels = false;\n');
end

if get(handles.ps_axis_at_origin ,'value') 
    fprintf(fid,'problem_space_display.axis_at_origin = true;\n');
else
    fprintf(fid,'problem_space_display.axis_at_origin = false;\n');
end

if get(handles.ps_axis_outside ,'value') 
    fprintf(fid,'problem_space_display.axis_outside_domain = true;\n');
else
    fprintf(fid,'problem_space_display.axis_outside_domain = false;\n');
end
if get(handles.ps_outer_boundaries ,'value') 
    fprintf(fid,'problem_space_display.outer_boundaries = true;\n');
else
    fprintf(fid,'problem_space_display.outer_boundaries = false;\n');
end
if get(handles.ps_cpml_boundaries ,'value') 
    fprintf(fid,'problem_space_display.cpml_boundaries = true;\n');
else
    fprintf(fid,'problem_space_display.cpml_boundaries = false;\n');
end
if get(handles.ps_xn_grid ,'value') 
    fprintf(fid,'problem_space_display.grid_xn = true;\n');
else
    fprintf(fid,'problem_space_display.grid_xn = false;\n');
end
if get(handles.ps_yn_grid ,'value') 
    fprintf(fid,'problem_space_display.grid_yn = true;\n');
else
    fprintf(fid,'problem_space_display.grid_yn = false;\n');
end
if get(handles.ps_zn_grid ,'value') 
    fprintf(fid,'problem_space_display.grid_zn = true;\n');
else
    fprintf(fid,'problem_space_display.grid_zn = false;\n');
end
if get(handles.ps_xp_grid ,'value') 
    fprintf(fid,'problem_space_display.grid_xp = true;\n');
else
    fprintf(fid,'problem_space_display.grid_xp = false;\n');
end
if get(handles.ps_yp_grid ,'value') 
    fprintf(fid,'problem_space_display.grid_yp = true;\n');
else
    fprintf(fid,'problem_space_display.grid_yp = false;\n');
end
if get(handles.ps_zp_grid ,'value') 
    fprintf(fid,'problem_space_display.grid_zp = true;\n');
else
    fprintf(fid,'problem_space_display.grid_zp = false;\n');
end

fclose(fid);

% ================================================================



function project_name_Callback(hObject, eventdata, handles)
% hObject    handle to project_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of project_name as text
%        str2double(get(hObject,'String')) returns contents of project_name as a double


% --- Executes during object creation, after setting all properties.
function project_name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to project_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in ps_toggle.
function ps_toggle_Callback(hObject, eventdata, handles)
col = [254 253 222]/255;
set(handles.ps_toggle,'value',0,'backgroundcolor',col);
set(handles.gsl_toggle,'value',1,'backgroundcolor','w');
set(handles.o_toggle,'value',1,'backgroundcolor','w');

set(handles.ps_panel,'visible','on');
set(handles.g_panel,'visible','off');
set(handles.sl_panel,'visible','off');
set(handles.o_panel,'visible','off');

% --- Executes on button press in gsl_toggle.
function gsl_toggle_Callback(hObject, eventdata, handles)
col = [254 253 222]/255;
set(handles.ps_toggle,'value',1,'backgroundcolor','w');
set(handles.gsl_toggle,'value',0,'backgroundcolor',col);
set(handles.o_toggle,'value',1,'backgroundcolor','w');

set(handles.ps_panel,'visible','off');
set(handles.g_panel,'visible','on');
set(handles.sl_panel,'visible','on');
set(handles.o_panel,'visible','off');

% --- Executes on button press in o_toggle.
function o_toggle_Callback(hObject, eventdata, handles)
col = [254 253 222]/255;
set(handles.ps_toggle,'value',1,'backgroundcolor','w');
set(handles.gsl_toggle,'value',1,'backgroundcolor','w');
set(handles.o_toggle,'value',0,'backgroundcolor',col);

set(handles.ps_panel,'visible','off');
set(handles.g_panel,'visible','off');
set(handles.sl_panel,'visible','off');
set(handles.o_panel,'visible','on');


% --- Executes on button press in TE_mode.
function TE_mode_Callback(hObject, eventdata, handles)
set(handles.TE_mode,'value',1);
set(handles.TM_mode,'value',0);
set(handles.TEM_mode,'value',0);
set(handles.E_phi,'enable','on');
set(handles.H_phi,'enable','off');
set(handles.E_x,'enable','off');
set(handles.E_y,'enable','off');
set(handles.kx_start,'enable','on');
set(handles.ky_start,'enable','on');

% --- Executes on button press in TM_mode.
function TM_mode_Callback(hObject, eventdata, handles)
set(handles.TE_mode,'value',0);
set(handles.TM_mode,'value',1);
set(handles.TEM_mode,'value',0);
set(handles.E_phi,'enable','off');
set(handles.H_phi,'enable','on');
set(handles.E_x,'enable','off');
set(handles.E_y,'enable','off');
set(handles.kx_start,'enable','on');
set(handles.ky_start,'enable','on');


function edit50_Callback(hObject, eventdata, handles)
% hObject    handle to edit50 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit50 as text
%        str2double(get(hObject,'String')) returns contents of edit50 as a double


% --- Executes during object creation, after setting all properties.
function edit50_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit50 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit51_Callback(hObject, eventdata, handles)
% hObject    handle to edit51 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit51 as text
%        str2double(get(hObject,'String')) returns contents of edit51 as a double


% --- Executes during object creation, after setting all properties.
function edit51_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit51 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ky_stop_Callback(hObject, eventdata, handles)
% hObject    handle to ky_stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ky_stop as text
%        str2double(get(hObject,'String')) returns contents of ky_stop as a double


% --- Executes during object creation, after setting all properties.
function ky_stop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ky_stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function kx_stop_Callback(hObject, eventdata, handles)
% hObject    handle to kx_stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of kx_stop as text
%        str2double(get(hObject,'String')) returns contents of kx_stop as a double


% --- Executes during object creation, after setting all properties.
function kx_stop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to kx_stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ky_step_Callback(hObject, eventdata, handles)
% hObject    handle to ky_step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ky_step as text
%        str2double(get(hObject,'String')) returns contents of ky_step as a double


% --- Executes during object creation, after setting all properties.
function ky_step_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ky_step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function kx_step_Callback(hObject, eventdata, handles)
% hObject    handle to kx_step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of kx_step as text
%        str2double(get(hObject,'String')) returns contents of kx_step as a double


% --- Executes during object creation, after setting all properties.
function kx_step_CreateFcn(hObject, eventdata, handles)
% hObject    handle to kx_step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in kxky.
function kxky_Callback(hObject, eventdata, handles)
set(handles.kxky,'value',1);
set(handles.kx_sweep,'value',0);
set(handles.ky_sweep,'value',0);
set(handles.kx_stop,'enable','off');
set(handles.kx_step,'enable','off');
set(handles.ky_stop,'enable','off');
set(handles.ky_step,'enable','off');

% --- Executes on button press in kx_sweep.
function kx_sweep_Callback(hObject, eventdata, handles)
set(handles.kxky,'value',0);
set(handles.kx_sweep,'value',1);
set(handles.ky_sweep,'value',0);
set(handles.kx_stop,'enable','on');
set(handles.kx_step,'enable','on');
set(handles.ky_stop,'enable','off');
set(handles.ky_step,'enable','off');

% --- Executes on button press in ky_sweep.
function ky_sweep_Callback(hObject, eventdata, handles)
set(handles.kxky,'value',0);
set(handles.kx_sweep,'value',0);
set(handles.ky_sweep,'value',1);
set(handles.kx_stop,'enable','off');
set(handles.kx_step,'enable','off');
set(handles.ky_stop,'enable','on');
set(handles.ky_step,'enable','on');

% --- Executes on button press in TM_mode.
function radiobutton6_Callback(hObject, eventdata, handles)
% hObject    handle to TM_mode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of TM_mode



function E_phi_Callback(hObject, eventdata, handles)
% hObject    handle to E_phi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of E_phi as text
%        str2double(get(hObject,'String')) returns contents of E_phi as a double


% --- Executes during object creation, after setting all properties.
function E_phi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to E_phi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function H_phi_Callback(hObject, eventdata, handles)
% hObject    handle to H_phi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of H_phi as text
%        str2double(get(hObject,'String')) returns contents of H_phi as a double


% --- Executes during object creation, after setting all properties.
function H_phi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to H_phi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in TE_mode.
function radiobutton7_Callback(hObject, eventdata, handles)
% hObject    handle to TE_mode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of TE_mode



function source_z_Callback(hObject, eventdata, handles)
% hObject    handle to source_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of source_z as text
%        str2double(get(hObject,'String')) returns contents of source_z as a double


% --- Executes during object creation, after setting all properties.
function source_z_CreateFcn(hObject, eventdata, handles)
% hObject    handle to source_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function reflection_z_Callback(hObject, eventdata, handles)
% hObject    handle to reflection_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of reflection_z as text
%        str2double(get(hObject,'String')) returns contents of reflection_z as a double


% --- Executes during object creation, after setting all properties.
function reflection_z_CreateFcn(hObject, eventdata, handles)
% hObject    handle to reflection_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in transmission_on.
function transmission_on_Callback(hObject, eventdata, handles)

if get(handles.transmission_on,'value')
    set(handles.transmission_z,'enable','on');
else
    set(handles.transmission_z,'enable','off');
end

function transmission_z_Callback(hObject, eventdata, handles)
% hObject    handle to transmission_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of transmission_z as text
%        str2double(get(hObject,'String')) returns contents of transmission_z as a double


% --- Executes during object creation, after setting all properties.
function transmission_z_CreateFcn(hObject, eventdata, handles)
% hObject    handle to transmission_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in reflection_on.
function reflection_on_Callback(hObject, eventdata, handles)
if get(handles.reflection_on,'value')
    set(handles.reflection_z,'enable','on');
else
    set(handles.reflection_z,'enable','off');
end


% --- Executes on button press in TEM_mode.
function TEM_mode_Callback(hObject, eventdata, handles)
set(handles.TE_mode,'value',0);
set(handles.TM_mode,'value',0);
set(handles.TEM_mode,'value',1);
set(handles.E_phi,'enable','off');
set(handles.H_phi,'enable','off');
set(handles.E_x,'enable','on');
set(handles.E_y,'enable','on');
set(handles.kx_start,'enable','off');
set(handles.ky_start,'enable','off');


function E_x_Callback(hObject, eventdata, handles)
% hObject    handle to E_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of E_x as text
%        str2double(get(hObject,'String')) returns contents of E_x as a double


% --- Executes during object creation, after setting all properties.
function E_x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to E_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function E_y_Callback(hObject, eventdata, handles)
% hObject    handle to E_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of E_y as text
%        str2double(get(hObject,'String')) returns contents of E_y as a double


% --- Executes during object creation, after setting all properties.
function E_y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to E_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
