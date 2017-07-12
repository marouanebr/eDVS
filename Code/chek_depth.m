function chek_depth(varargin)
% CHECK_DEPTH Return the 3D real world coordinates of an event in cm unit
%
%   Usage  
%   chek_depth   
%   Fetch the necessary input data files from the folder /take
%
%   chek_depth(tol)
%   Specify the tolerance value i.e. the radius of the edges
%
%   chek_depth(tol,'path')
%   Fetch the necessary input data files from the folder /path
%
%   Description
%   The function fetchs the following data from the necessary files
%   - Coordinates of the markers on the eDVS camera
%   - Coordinates of the markers on the top vertices of the boxes
%   - The events time stamps and positions on the sensor chip
%   Results:
%   - A figure with the take reconstraction 
%   - A n*3 matrix result, where n is the number of events:
%     result: event time stamp | computed depth | real depth
%   - A n*7 matrix result_viz, used in script to for visualization
%     result_viz is just an extension of the matrix result   
%     result_viz: result | event's coordinates (xyz) | tracking time stamp
%  
%   Example
%   chek_depth(5,'C:\Users\User_Name\Downloads\take')
%
%   The unit used for all variables is cm
%   See also
%   draw_boxes(), correct_camera(), events_viz()
% ------
% Author: Marouane Ben Romdhane
% e-mail: marouanebr@gmail.com
% Created: 01.03.2016,    using Matlab R2015b


%% configure the parameters of the figure 
figure;
title('scene reconstruction');
grid on; hold on;
xlabel('x'); ylabel('y'); zlabel('z'); 
%axis([-200 200 -200 200 -200 200]);

set(gcf, 'renderer', 'opengl');
%newpos = [-50 100 0]; target = [-50 50 0];
%set(gca, 'CameraPosition', newpos, 'CameraTarget', target);

%% initialization (check REPORT.PDF for more details) 
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
% load and draw the boxes
[n,top,bottom] = draw_boxes(tol,path);

% fetch the camera data in matrix S
camera_file = [path{:},'\Camera.log'];  
S = importdata(camera_file,' ');

% fetch the event data in matrix T
event_file = [path{:},'\dvs_events.tsv'];
T = dlmread(event_file,'\t',3, 0);
T(size(T,1)+1,:) = -1; % add a row of -1 at the end of T
ts = find(T(:,2)== -1); % create a matrix with the indices of all the 
                        % tracking software time stamps
[q,~] = size(ts);

% the matrix holding the final result
result = zeros(0,3); % event time stamp | computed depth | real depth
result_viz = zeros(0,7); % used to visualize the results 

%% compute results for each new time stamp of the tracking system
for i = 1:q-1 
    % ignore the cases when there are no new events  
    if ts(i) ~= (ts(i+1)-1)        
        %% load a new event matrix
        % get the number of events for the new time stamp of the tracking system
        [s,~] = size(T((ts(i)+1):(ts(i+1)-1),2:3)); 
        event = zeros(s,4); % event = px*py*z*id
        event(:,1:2) = T((ts(i)+1):(ts(i+1)-1),2:3); % px and py 
        event(:,3) = zeros(s,1); % missing depth
        event(:,4) = T((ts(i)+1):(ts(i+1)-1),1); % event time stamp    

        %% load a new camera matrix
        % check for the correct time stamp of the tracking system
        [~, index] = min(abs(S(:,1)-  T(ts(i),1)));          
        cam = zeros(3,3);
        % load the camera coordinates
        % *100 to change to cm
        for j=1:3    
            cam(j,:)= 100 .* S(index,17+6*j:19+6*j);   
        end

        %% correct the camera position
        cam = correct_camera(cam);
        drawPoint3d(cam); % draw the camera  
        drawPolygon3d(cam,'k');        

        %% compute the result matrix
        [a,px,py] = field_of_view(cam,los,fov); % get the field of view 
        % create a matrix with the projections of all events on the FOV
        b = repmat(a,s,1);        
        proj = b + event(:,1)*px + event(:,2)*py;

        %% check for intersections between the lines and the boxes
        % initialization for each new event matrix
        % res matrix: event time stamp|depth|real depth|xyz|tracking ts
        res = zeros(s,7);
        % event time stamp
        res(:,1) = event(:,4);     
        % event depth to be tested
        res(:,2) = event(:,3);
        % time stamp of the tracking system
        res(:,7) = S(index,1);
        % D holds the depth of each intersection point        
        D = zeros(n,12); % D = box * depth from each edge
        % PT holds the intesections coordinates        
        PT = zeros(n,12,3); % PT = box * edge * coordinates
        % edge holds edges of the box
        edge = zeros(12,6); % edge = edge number * [top_point, bottom_point]

        % check for intersection between the lines and the boxes
        % and return the intersection point
        for j=1:s % for each event
            % create a line crossing the camera and the event projection
            L = createLine3d(cam(1,:), proj(j,:)); 
            for k=1:n % for each box      
                % create a matrix CYL holding all the box's edges as cylinder    
                for l=1:4            
                    % vertical edges
                    edge(l,:) = [permute(top(k,l,:),[1 3 2]) permute(bottom(k,l,:),[1 3 2])];
                    % top horizontal edges
                    edge(l+4,:) = [permute(top(k,l,:),[1 3 2]) permute(top(k,l+1,:),[1 3 2])];          
                    %bottom horizontal edges
                    edge(l+8,:) = [permute(bottom(k,l,:),[1 3 2]) permute(bottom(k,l+1,:),[1 3 2])];
                end                       
                CYL = [edge repmat(tol,12,1)];
                % intersection with each edge
                for temp = 1:12
                    inter = intersectLineCylinder(L, CYL(temp,:));                    
                    if size(inter,1) == 0; % no intersection
                        inter = [NaN, NaN, NaN];
                    end 
                    PT(k,temp,:) = inter(1,:);
                    point = permute(PT(k,temp,:),[1 3 2]);   
                    D(k,temp) = distancePoints3d(cam(1,:), point);                     
                end  
            end
            
            % check for the closest, non zero, intersection point 
            [min_val,indice] = nanmin(D(:)); 
            [row,col] = ind2sub(size(D),indice);

            % draw the intersection coordinates and save the real depth
            if(~isnan(min_val))
                drawPoint3d(permute(PT(row,col,:),[1 3 2]));             
                res(j,3)= min_val; % save the depth in 3rd column
                res(j,4:6) = permute(PT(row,col,:),[1 3 2]); % save the xyz
            else
                res(j,3)= NaN;
                res(j,4:6) =[NaN, NaN, NaN];
            end   
        end
        result = vertcat(result,res(:,1:3));        
        result_viz = vertcat(result_viz,res);
        disp('***************  new time stamp of the tracking system ***************'); 
    end    
end
% save the final matrix
save('result','result');
save('result_viz','result_viz');

% save the figure
print(title,'-dpng');
savefig(title);
end
