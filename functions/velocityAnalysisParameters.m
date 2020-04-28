
function[out] = velocityAnalysisParameters
function [varargout] = cmpGUI(varargin)
% CMPGUI MATLAB code for cmpGUI.fig
%      CMPGUI, by itself, creates a new CMPGUI or raises the existing
%      singleton*.
%
%      H = CMPGUI returns the handle to a new CMPGUI or the handle to
%      the existing singleton*.
%
%      CMPGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CMPGUI.M with the given input arguments.
%
%      CMPGUI('Property','Value',...) creates a new CMPGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before cmpGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to cmpGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help cmpGUI

% Last Modified by GUIDE v2.5 27-Apr-2020 22:18:51

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @cmpGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @cmpGUI_OutputFcn, ...
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

end

% --- Executes just before cmpGUI is made visible.
function cmpGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to cmpGUI (see VARARGIN)

% Choose default command line output for cmpGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes cmpGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);
end

% --- Outputs from this function are returned to the command line.
function varargout = cmpGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
end


function ndir_Callback(hObject, eventdata, handles)
% hObject    handle to ndir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ndir as text
%        str2double(get(hObject,'String')) returns contents of ndir as a double
end

% --- Executes during object creation, after setting all properties.
function ndir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ndir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function nref_Callback(hObject, eventdata, handles)
% hObject    handle to nref (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nref as text
%        str2double(get(hObject,'String')) returns contents of nref as a double
end

% --- Executes during object creation, after setting all properties.
function nref_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nref (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
closereq();
return
end


% --- Executes when selected object is changed in uibuttongroup3.
function uibuttongroup3_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uibuttongroup3 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
L1 = handles.L1.Value;
L2 = handles.L2.Value;
disp('norm')
end
end