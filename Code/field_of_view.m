function [a,px,py] = field_of_view(cam,los,fov)

% field_of_view Draws the field of view 
%
%   Usage  
%   [a, px, py] = field_of_view(cam)   
%   Fetch the coordinates of camera's tracking balls
%   
%   Description
%   The function fetchs the coordinates of the tracking balls on the eDVS.
%   Then it computes the coordinates of the corners of the Field of View. 
%   Finally, it returns the following variables:
%   1- A: 1*3 vector containing the coordinates of the point A
%   2- px: 1*3 vector corresponds to the vector AB/128, 128 corresponds to
%          resolution of the eDVS
%   3- py: 1*3 vector corresponds to the vector AD/128, 128 corresponds to
%          resolution of the eDVS
%
%                         A---------B
%                         |         |
%                         E    H    |
%                         |         |
%                         D---------C
%   
%
%   See also
%
% ------
% Author: Marouane Ben Romdhane
% e-mail: marouanebr@gmail.com
% Created: 05.12.2015,    using Matlab R2015b

%% compute 'h', the projection of the eDVS on the FOV
% compute the projection of the back tracking ball of the eDVS, on the line
% created by the other two tracking balls on the eDVS
line = createLine3d(cam(1,:),cam(2,:));
proj_1 = projPointOnLine3d(cam(3,:), line);

% vector N colinear to the Line of Sight, i.e. perpendicular to the FOV
N = proj_1 - cam(3,:);
% normalization of the vector N:= norm(n)=los
N = (los/norm(N))*N; 

% the projection of camera on the FOV i.e. the center of the FOV
h = cam(1,:)+N;

%% find the coordinates of the corners of the FOV
% normalized vector u corresponds to vector AB, i.e. the side of the FOV
u = (cam(2,:)-cam(1,:));
u = (fov/norm(u))*u;

% compute the coordinates of E
e = h-(u/2);

% normalized vector v corresponds to vector AD
v = -cross(N,u);
v = (fov/norm(v))*v;

a = e+(v/2); 
b = a+u; 
c = b-v; 
d = a-v; 

%% return px and py
% vectors px and py corresponds to a horizontal/vertical one pixel step
% on the camera's lens
px = u/128;
py = -v/128;

%% draw the FOV
% compute and draw the 4 corners of the FOV

%draw the point h and the edge [cam-h]
drawEdge3d([cam(1,:) h]);

% draw the FOV's edges
drawPoint3d(a); drawPoint3d(b); drawPoint3d(c); drawPoint3d(d);
drawEdge3d([a b]); drawEdge3d([b c]); 
drawEdge3d([c d]); drawEdge3d([a d]); 

