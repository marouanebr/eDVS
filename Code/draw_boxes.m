function [n,top,bottom] = draw_boxes(tol,path)
% DRAW_BOXES Draws the boxes from the test enviorment 
%
%   Usage
%   [n, top, bottom] = DRAW_BOXES(tol,path)
%
%   Description
%   tol: the radius of a cylinder around each edge
%   path: Location of the folder containing the log files, The log file
%   name of tracking balls of each box should have this format:
%   'Box(i).log' where (i) is in [1,..,n] %            
%
%   the dimension of the box matrices is n*m*p where:
%       - n : the box's index
%       - m : the vertex's index
%       - p : the coordinate's index (x, y or z)
%
%
%   The function do the following steps:
%   1- fetch the coordinates of the top vertices of each Box i and save 
%      them in the matrix top()
%   2- rearrange the vertices to be clockwise (CW) or counterclockwise (CCW) 
%   3- create a mirror matrix for the bottom vertices bottom() 
%   4- draw the boxes 
%   5- draws a cylinder with radius = tolerance around each edge
%   6- return the matrix top and bottom
% 
%   The unit used for all variables is cm
%   See also
%   chek_depth(), correct_camera(), events_viz()
% ------
% Author: Marouane Ben Romdhane
% e-mail: marouanebr@gmail.com
% Created: 01.03.2016,    using Matlab R2015b

%% determine the number of the boxes i.e. log files
n = 1;
str = [path{:},'\Box',num2str(n),'.log'];

while(exist(str,'file')==2)    
    n=n+1;
    str = ['Box',num2str(n),'.log'];
end    

n=n-1;    
str_1=['number of the boxes is ',num2str(n)];
disp(str_1);

%% create the top matrix
top=zeros(n,4,3);

for i=1:n
    % fetch from the log file the mean value of the first 100 lines
    filename=[path{:},'\Box',num2str(i),'.log'];
    delimiterIn = ' '; 
    % *100 to change the unit in cm
    A = 100 .* importdata(filename,delimiterIn);
    for j=1:4
        top(i,j,:)= mean(A(1:100,17+6*j:19+6*j)); 
    end      
     
    %% Rerrange the marks CW or CCW
    %2D transformation of the top face
    x = top(i,:,1);
    y = top(i,:,2);
    z = top(i,:,3);
    cx = mean(x);
    cz = mean(z);
    
    % sort the vertices of the top CW or CCW
    a = atan2(z-cz, x-cx);
    [~,order] = sort(a);
    top(i,:,1) = x(order);
    top(i,:,2) = y(order);
    top(i,:,3) = z(order);          
    
end

% add a 5th vertice, replicating the first one (useful to draw the boxes)
top(:,5,:) = top(:,1,:);

% create a matrix similar to the matrix "top" for the bottom vertices
bottom = top;
bottom(:,:,2) = 0;                   

%% draw the boxes
% draw the side faces
for i=1:n
    % draw the side faces
    for j=1:4
        POLY = [permute(top(i,j,:),[1 3 2]);...
                  permute(top(i,j+1,:),[1 3 2]);...
                  permute(bottom(i,j+1,:),[1 3 2]);...
                  permute(bottom(i,j,:),[1 3 2])];
        fillPolygon3d(POLY,'y');
    end
    
    % draw the top face
	POLY=[permute(top(i,1,:),[1 3 2]);...
            permute(top(i,2,:),[1 3 2]);...
            permute(top(i,3,:),[1 3 2]);...
            permute(top(i,4,:),[1 3 2])];
    fillPolygon3d(POLY,'y');    
end

% draw a gray ground plane
plane = createPlane(permute(bottom(1,1,:),[1 3 2]),...
                    permute(bottom(1,2,:),[1 3 2]),...
                    permute(bottom(1,3,:),[1 3 2]));
drawPlane3d(plane,[0.7 0.7 0.7])

%% draw cylinders
edge = zeros(12,6); % edge = edge number * [top_point, bottom_point]

for i=1:n % for each box      
    % create a matrix CYL holding all the box's edges as cylinder    
    for j=1:4            
        % vertical edges
        edge(j,:) = [permute(top(i,j,:),[1 3 2]) permute(bottom(i,j,:),[1 3 2])];
        % top horizontal edges
        edge(j+4,:) = [permute(top(i,j,:),[1 3 2]) permute(top(i,j+1,:),[1 3 2])];          
        %bottom horizontal edges
        edge(j+8,:) = [permute(bottom(i,j,:),[1 3 2]) permute(bottom(i,j+1,:),[1 3 2])];
    end                       
    CYL = [edge repmat(tol,12,1)];
    for j = 1:12
        drawCylinder(CYL(j,:));                     
    end  
end

