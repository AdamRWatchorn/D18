% 3D Rotation around the Y axis
%
% x is a matrix with points to rotate th is the rotation angle in radians (positive is CCW)
%
% The function returns the transformed points, and the transformation matrix

function [x_ry, Mry]=Ry(x,th)

 % 3D homogeneous translation matrix
 Mry=[cos(th)   0  sin(th) 0
       0        1     0    0
      -sin(th)  0  cos(th) 0
       0        0     0    1];
     
 x_ry=Mry*x;
