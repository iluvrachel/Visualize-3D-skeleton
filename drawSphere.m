function [Sphere] = drawSphere(pos,r)
%DRAWSPHERE 此处显示有关此函数的摘要
%   此处显示详细说明
x0 = pos(1);
y0 = pos(2);
z0 = pos(3);
[x,y,z]=sphere(19);
Sphere = mesh(x0+r*(-x),y0+r*(-y),z0+r*(-z));
set(Sphere,'FaceColor',[0,0,0]/255);
set(Sphere,'EdgeAlpha',0);
% Sphere.FaceNormalsMode = 'manual';
% Sphere.FaceNormals = -Sphere.FaceNormals; 
end

