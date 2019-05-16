function VisSkeleton_obj(position)

position=position(:,:,1:4:end);
frameCount=size(position,3);

points_1_id=[2,11,12,13,14]; %body
points_2_id=[18,17,16,15,12,19,20,21,22]; %upper limbs
points_3_id=[6,5,4,3,2,7,8,9,10]; %lower limbs

points_id={points_1_id,points_2_id,points_3_id}';

joints_to_show=[2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22];

x=position(1,:,:);
y=position(2,:,:);
z=position(3,:,:);
max_point=[max(x(:)),max(y(:)),max(z(:))];
min_point=[min(x(:)),min(y(:)),min(z(:))];

mesh_X = [];
mesh_Y = [];
mesh_Z = [];

node_mesh_X = [];
node_mesh_Y = [];
node_mesh_Z = [];


for i=1:frameCount
    disp(i-1);
    skeleton_p=position(:,:,i)';  %(23,3)
    fig=figure();
    hold on
    
    plot3(-2,-2,-0.5,'.white');
    plot3(-2, 2,-0.5,'.white');
    plot3(2,-2,-0.5,'.white');
    plot3(2,2,-0.5,'.white'); 
    plot3(-2,-2,2.5,'.white');
    plot3(-2,2,2.5,'.white');
    plot3(2,-2,2.5,'.white');
    plot3(2,2,2.5,'.white');
    view(90,5);  %Figure1
    axis equal
    axis tight
    axis off
    
    for j=1:size(points_id,1)
        p_id=points_id{j};
        points=skeleton_p(p_id,:);
        
        for p=2:size(points,1)
            cy = drawCylinder(points(p-1,:)',points(p,:)',0.2,10,[168,37,37]/255,1,0);
            % get xyz data
            mesh_X = cat(3,mesh_X,cy.XData);
            mesh_Y = cat(3,mesh_Y,cy.YData);
            mesh_Z = cat(3,mesh_Z,cy.ZData);

        end
    end
    
    for p=1:size(joints_to_show,2)
        cur_p=joints_to_show(p);
        points=skeleton_p(cur_p,:);
        % magnify head(joint number 14)
        if cur_p==14
            sp = drawSphere(points(:,:)',0.6);
        else
            sp = drawSphere(points(:,:)',0.3);
        end
        
        node_mesh_X = cat(3,node_mesh_X,sp.XData);
        node_mesh_Y = cat(3,node_mesh_Y,sp.YData);
        node_mesh_Z = cat(3,node_mesh_Z,sp.ZData);

    end

    
    hold off

    fid=fopen('test.obj','w');
    nn=1;
    [nCy,mm] = WriteV(fid,mesh_X,mesh_Y,mesh_Z,nn);
    [nSp,useless] = WriteV(fid,node_mesh_X,node_mesh_Y,node_mesh_Z,mm);
    fprintf(fid,'g mesh\n');
    WriteF(fid,mesh_X,mesh_Y,mesh_Z,nCy);
    WriteF(fid,node_mesh_X,node_mesh_Y,node_mesh_Z,nSp);
    fprintf(fid,'g\n\n');
    fclose(fid);

    
    close all

end

end


% write v $ vt
% fid is the .obj file to write
% x y z is position of each vertice
% nn is the current order of vertice
function [n,m]= WriteV(fid,x,y,z,nn)
normals=0;
l=size(x,1); h=size(x,2); d=size(x,3);

n=zeros(l,h,d);

for k=1:d
    for i=1:l
        for j=1:h
            n(i,j,k)=nn;
            fprintf(fid, 'v %f %f %f\n',x(i,j,k),y(i,j,k),z(i,j,k));
            fprintf(fid, 'vt %f %f\n',(i-1)/(l-1),(j-1)/(h-1));
            if (normals) fprintf(fid, 'vn %f %f %f\n', nx(i,j,k),ny(i,j,k),nz(i,j,k)); end
            nn=nn+1;
        end
    end
end
m = nn;
end

% write f 
% fid is the .obj file to write
% x y z is position of each vertice
% n is a 2D array return from WriteV
function WriteF(fid,x,y,z,n)
normals=0;%(ignore normals)
l=size(x,1); h=size(x,2); d=size(x,3);
for k=1:d
    for i=1:(l-1)
        for j=1:(h-1)
            if (normals)
                fprintf(fid,'f %d/%d/%d %d/%d/%d %d/%d/%d %d/%d/%d\n',n(i,j,k),n(i,j,k),n(i,j,k),n(i+1,j,k),n(i+1,j,k),n(i+1,j,k),n(i+1,j+1,k),n(i+1,j+1,k),n(i+1,j+1,k),n(i,j+1,k),n(i,j+1,k),n(i,j+1,k));
            else
                fprintf(fid,'f %d/%d %d/%d %d/%d %d/%d\n',n(i,j,k),n(i,j,k),n(i+1,j,k),n(i+1,j,k),n(i+1,j+1,k),n(i+1,j+1,k),n(i,j+1,k),n(i,j+1,k));
            end
        end
    end
end

end


function [Sphere] = drawSphere(pos,r)
%   x0,y0,z0 is the centre of sphere
x0 = pos(1);
y0 = pos(2);
z0 = pos(3);
[x,y,z]=sphere(19);
% here use minus value to keep the normals of each Face outward
Sphere = mesh(x0+r*(-x),y0+r*(-y),z0+r*(-z));
% black
set(Sphere,'FaceColor',[0,0,0]/255);
% disable edges
set(Sphere,'EdgeAlpha',0);
% Sphere.FaceNormalsMode = 'manual';
% Sphere.FaceNormals = -Sphere.FaceNormals; 
end



function [Cylinder EndPlate1 EndPlate2] = drawCylinder(X1,X2,r,n,cyl_color,closed,lines)
%
% This function constructs a cylinder connecting two center points 
% 
% Usage :
% [Cylinder EndPlate1 EndPlate2] = drawCylinder(X1+20,X2,r,n,'r',closed,lines)
%    
% Cylinder-------Handle of the cylinder
% EndPlate1------Handle of the Starting End plate
% EndPlate2------Handle of the Ending End plate
% X1 and X2 are the 3x1 vectors of the two points
% r is the radius of the cylinder
% n is the no. of elements on the cylinder circumference (more--> refined)
% cyl_color is the color definition like 'r','b',[0.52 0.52 0.52]
% closed=1 for closed cylinder or 0 for hollow open cylinder
% lines=1 for displaying the line segments on the cylinder 0 for only
% surface
% 
% Typical Inputs
% X1=[10 10 10];
% X2=[35 20 40];
% r=1;
% n=20;
% cyl_color='b';
% closed=1;
%

%see more information please go to www.matlabsky.cn


% NOTE: There is a MATLAB function "cylinder" to revolve a curve about an
% axis. This "Cylinder" provides more customization like direction and etc


% Calculating the length of the cylinder
length_cyl=norm(X2-X1);

% Creating a circle in the YZ plane
t=linspace(0,2*pi,n)';
x2=r*cos(t);
x3=r*sin(t);

% Creating the points in the X-Direction
x1=[0 length_cyl];

% Creating (Extruding) the cylinder points in the X-Directions
xx1=repmat(x1,length(x2),1);
xx2=repmat(x2,1,2);
xx3=repmat(x3,1,2);

% Drawing two filled cirlces to close the cylinder
if closed==1
    hold on
    EndPlate1=fill3(xx1(:,1),xx2(:,1),xx3(:,1),'r');
    EndPlate2=fill3(xx1(:,2),xx2(:,2),xx3(:,2),'r');
end

% Plotting the cylinder along the X-Direction with required length starting
% from Origin
Cylinder=mesh(xx1,xx2,xx3);

% Defining Unit vector along the X-direction
unit_Vx=[1 0 0];

% Calulating the angle between the x direction and the required direction
% of cylinder through dot product
angle_X1X2=acos( dot( unit_Vx,(X2-X1) )/( norm(unit_Vx)*norm(X2-X1)) )*180/pi;

% Finding the axis of rotation (single rotation) to roate the cylinder in
% X-direction to the required arbitrary direction through cross product
axis_rot=cross([1 0 0],(X2-X1) );

% Rotating the plotted cylinder and the end plate circles to the required
% angles
if angle_X1X2~=0 % Rotation is not needed if required direction is along X
    rotate(Cylinder,axis_rot,angle_X1X2,[0 0 0])
    if closed==1
        rotate(EndPlate1,axis_rot,angle_X1X2,[0 0 0])
        rotate(EndPlate2,axis_rot,angle_X1X2,[0 0 0])
    end
end

% Till now cylinder has only been aligned with the required direction, but
% position starts from the origin. so it will now be shifted to the right
% position
if closed==1
    set(EndPlate1,'XData',get(EndPlate1,'XData')+X1(1))
    set(EndPlate1,'YData',get(EndPlate1,'YData')+X1(2))
    set(EndPlate1,'ZData',get(EndPlate1,'ZData')+X1(3))
    
    set(EndPlate2,'XData',get(EndPlate2,'XData')+X1(1))
    set(EndPlate2,'YData',get(EndPlate2,'YData')+X1(2))
    set(EndPlate2,'ZData',get(EndPlate2,'ZData')+X1(3))
end
set(Cylinder,'XData',get(Cylinder,'XData')+X1(1))
set(Cylinder,'YData',get(Cylinder,'YData')+X1(2))
set(Cylinder,'ZData',get(Cylinder,'ZData')+X1(3))

% Setting the color to the cylinder and the end plates
set(Cylinder,'FaceColor',cyl_color)
if closed==1
    set([EndPlate1 EndPlate2],'FaceColor',cyl_color)
else
    EndPlate1=[];
    EndPlate2=[];
end

% If lines are not needed making it disapear
if lines==0
    set(Cylinder,'EdgeAlpha',0)
end



