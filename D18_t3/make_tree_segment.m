
function []=make_tree_segment(shape,Mp,phi,depth)
 
 % A 2D box centered along the Y axis
 %x=[-.5:.01:.5]';
 %y=[0:.01:1]';
 %box=[x                   zeros(length(x),1)
 %    -.5*ones(length(x),1)       y
 %     x                   ones(length(x),1)
 %     .5*ones(length(x),1)        y           ];
 %box(:,end+1)=zeros(length(box),1);
 %box(:,end+1)=ones(length(box),1);
 
 % We will use rectangular parts, we know we don't want to
 % use non-uniform scaling with hierarchical transforms, so
 % we will prescale the box
 bT=S(shape',[1 1 1 1]');% bT is not a rectangular 'part'
 % that is 2 units tall and 1 unit
 % wide
 
 
 [temp,Ms]=S([0 0 0 0]',[.8 .8 1 1]');
 % Now apply any rotation needed. In this case, rotation in the x-y plane by
 % angle theta means a rotation around the z axis by theta.
 [temp,Mr]=Rz([0 0 0 0]',phi);
 
 % Finally, translate to the top of part1, since part1 is 2 units long
 % along y, simply translate by 2 along y
 [temp,Mt]=T([0 0 0 0]',[0 1 0 1]');
 
 ML=eye(4,4);% The identity matrix to represent the first part's transform
 
 ML = Mt *Mr * Ms;
 % Part 1 is just M1*bT so
 b=Mp*ML*bT;
 
% figure(1);clf;hold on;grid on;axis equal;
% title('Hierarchical transforms example');
 plot(b(1,:),b(2,:),'b.');
 
 phi2 = 45;
 depth=depth+1;
 if(depth < 10)
 make_tree_segment(shape,ML,phi2,depth);
 end
end