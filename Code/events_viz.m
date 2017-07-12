function events_viz(varargin)
% EVENTS_VIZ Return the 3D real world coordinates of an event in cm unit
%
%   Usage  
%   events_viz   
%   Fetch the necessary input data files from the folder /take
%
%   events_viz(tol)
%   Specify the tolerance value
%
%   events_viz(tol,'path')
%   Fetch the necessary input data files from the folder /path
%
%   Description
%   The function fetchs the following data from the necessary files
%   - Coordinates of the tracking balls on the eDVS camera
%   - Coordinates of the top vertices of the boxes
%   - Data from the eDVS of the event
%   The user gives the time stamp of the first event to check, followed by
%   the number of the events. The function displays a figure with:
%   - the boxes in the scene
%   - the eDVS at the corresponding time
%   - the field of view
%   - the projection of the events on the edge of the boxes
%   - a line crossing the left lens and the event projection
%   Finally, a .png picture with the name "Event Visualization" is saved
%   
%   Example
%   events_viz('C:\Users\User_Name\Downloads\take')
%
%   The unit used for all variables is cm
%   See also
%   draw_boxes(), correct_camera(), chek_depth()
% ------
% Author: Marouane Ben Romdhane
% e-mail: marouanebr@gmail.com
% Created: 01.03.2016,    using Matlab R2015b

%% initialization (check README.PDF for more details)
clc; 
los = 26.067; % line of sight of the eDVS
fov = 21; % side of the Field of View (FOV) of the eDVS
% specify the tolerance and the folder for the necessary files
tol = 3; % corresponds to the radius of the edges 
path = {'take'};
if nargin~=0
    tol = varargin{1};
    if length(varargin) == 2
        path = varargin(2);
    end
end

%% fetch data from necessary files
% fetch the camera data in matrix S
camera_file = [path{:},'\Camera.log'];
S = importdata(camera_file,' ');
cam = zeros (3,3);
%test_mat: event_id | computed depth | real depth | xyz | eDVS time_stamp
load('result_viz','result_viz');

%% start the process of drawing the events
rep = 'y'; % confirms checking new series of events to test
rep_numb = 0; % number or repetitions

while (strcmp(rep,'y')||strcmp(rep,'yes'))
    %% configure the parameters of the figures
    figure;
    rep_numb = rep_numb+1;
    figure_name = ['Events Visualization ', num2str(rep_numb)];
    title(figure_name);
    grid on; hold on;
    xlabel('x'); ylabel('y'); zlabel('z'); 
    %axis([-150 200 -50 150 -200 150]);
    set(gcf, 'renderer', 'opengl');

    %% draw the boxes
    draw_boxes(tol,path);
   
    %% events to check
    % read event id
    str = 'Insert the time stamp of the first event ';
    id = input(str);
    
    % read the number of events
    str_1 = 'Insert the number of events to be checked ';
    i = input(str_1);
    
    %% draw the camera and the FOV
    % find the time stamp of the tracking system
    % find the index of the event id in the result_viz matrix
    [id_viz,~] = find(id == result_viz(:,1));
    % load the camera position at the correspondant time stamp
    [id_ts,~] = find(S(:,1) == result_viz(id_viz,7));      
    for j=1:3    
        cam(j,:)=100*S(id_ts,17+6*j:19+6*j);   
    end
    % correct the camera position
    cam = correct_camera(cam);
    drawPoint3d(cam); % draw the camera  
    drawPolygon3d(cam,'k');
    
    % draw the field of view
    field_of_view(cam,los,fov); % get the field of view 
   
    %% draw only the events on the edges
    % fill the event matrix: xyz of each event
    event = result_viz(id_viz:id_viz+i-1,4:6);
    % draw events and lines
    for j = 1:i
        if ~isnan(event(j,1))
            drawPoint3d(event(j,:));
            L = createLine3d(cam(1,:), event(j,:)); 
            drawLine3d(L);
        end
    end
    % save the figure
    print(figure_name,'-dpng');
    str_2 = 'do yo want to repeat ? (write ''y'' or ''yes'' to confirm )  ';
    rep = input(str_2,'s');
end