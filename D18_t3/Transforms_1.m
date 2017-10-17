% Simple example of transformations in 2D and 3D

% Let's start with a 2D box
x=[0:.01:1]';
y=[0:.01:1]';

% This makes a box with corners (0,0) - (0,1) - (1,1) - (1,0)
% Note that since we have 3D transform functions, we're going
% to take advantage of the fact that 2D points are just 3D
% points with z=0. So in the array below:
% The first column contains x values
% The second column contains y values
% The third column is all zeros (corresponding to z=0 for all points)
% and the last column is all ones (corresponding to the homogeneous coordinate w=1)
box=[x                   zeros(length(x),1)
     zeros(length(x),1)       y
     x                   ones(length(x),1)
     ones(length(x),1)        y           ];
box(:,end+1)=zeros(length(box),1);
box(:,end+1)=ones(length(box),1);
     
figure(1);clf;
title('A simple box');axis equal; hold on; grid on;
plot(box(:,1),box(:,2),'b.');

% Now, let's transform this box a bit.
% First example - Scale the box by (2,1) (make a chubby rectangle), 
%  rotate by pi/4 radians, and translate to (3,2)

bT=S(box',[2 1 1 1]');  % Note that transforms require points to be
			% column vectors, so we have to pass box'.
			% note the scale vector in homogeneous
			% coordinates.
			
% bT has a scaled box at this point, rotate it. Rotating in the XY
% plane is just the same as rotating around the Z axis in 3D
bT=Rz(bT,pi/4);

% Not bT has a scaled, rotated box. Translate it to the desired
% spot
bT=T(bT,[3 2 0 1]');

% Plot the transformed box in a different color (red for this example)
plot(bT(1,:), bT(2,:), 'r.');


% Second example - showing that the order of transformations matters!
% Same transformstions as above, but apply: Scale, Translate, Rotate
bT=S(box',[2 1 1 1]');
bT=T(bT,[3 2 0 1]');
bT=Rz(bT,pi/4);

% Plot this in magenta
plot(bT(1,:), bT(2,:), 'm.');

% See if you can implement the following:
% - From the starting box, set up a sequence of transformations that
%   transforms the box into an elongated diamond shape (twice as 
%   wide as it is tall), and located 4 units away from the origin
%   along a line 135 degrees away from horizontal.
% Plot it in any color you like.


%bT=S(box',[1 2 1 1]');

%bT=Rz(box',pi/4);

%bT=S(bT,[2 1 1 1]');

%bT=T(bT,[3 2 0 1]');

%plot(bT(1,:), bT(2,:), 'g.');

bT=T(box', [-.5 -.5 0 1]');
bT=Rz(bT, pi/4);
bT=S(bT, [2 1 1 1]');
bT=T(bT, [4 0 0 1]');
bT=Rz(bT, (135/360)*2*pi);

plot(bT(1,:), bT(2,:), 'g.');