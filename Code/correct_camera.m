function cam = correct_camera(camera)
% CORRECT_CAMERA Corrects the camera position
%
%   Usage
%   cam = CORRECT_CAMERA(camera)
%
%   the dimension of the camera matrix is n*m = 3*3 where:
%       - n : the index of the marker on the camera
%       - m : the coordinate's index (x, y or z)
%
%   The function do the following steps:
%   1- fetch the coordinates of the markers of the camera 
%   2- rearrange the vertices to be in the following order
%         I- master (left) sensor
%         II - slave (right)sensor 
%         III- back tracking ball
%      the user should put the makers exactly above the master 
%      sensor (check manual for more details) 
%   3- correct the height of the three markers
%   4- return the matrix cam
% 
%   The unit used for all variables is cm
%   See also
%   chek_depth(), events_viz()
% ------
% Author: Marouane Ben Romdhane
% e-mail: marouanebr@gmail.com
% Created: 01.03.2016,    using Matlab R2015b

%% intialization
cam=zeros(3,3);

%% coorect the tracking balls order
% create a matrix D holding the sum of length of the two edges next to each
% vertex
D = zeros(3);
D(1) = distancePoints3d(camera(1,:),camera(2,:)) + ...
            distancePoints3d(camera(1,:),camera(3,:));
D(2) = distancePoints3d(camera(2,:),camera(1,:)) + ...
            distancePoints3d(camera(2,:),camera(3,:));
D(3) = distancePoints3d(camera(3,:),camera(1,:)) + ...
            distancePoints3d(camera(3,:),camera(2,:));


% sort the vertex depending on their relative value in D
% master lens has the greatest value
% slave lens comes second
% back tracking ball comes last
[~,I] = sort(D,'descend');
cam(1,:) = camera(I(1),:);
cam(2,:) = camera(I(2),:);
cam(3,:) = camera(I(3),:);     

%% correct the height of the tracking balls
u = cam(2,:) - cam(1,:);
v = cam(3,:) - cam(2,:);
w = cross(u,v);
w = (1.5/norm(w))*w; % set the vertical shift value here
cam(1,:) = cam(1,:) + w;
cam(2,:) = cam(2,:) + w;
cam(3,:) = cam(3,:) + w;


