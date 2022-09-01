% test rotation

x0 = 1; % based on horizontal and vertical cartesian system
y0 = sqrt(3);
figure
axis equal
xran = [-4 4];
yran = [-4 4];
plot([0 0],[-4 4],'k--'); hold on
plot([-4 4],[0 0],'k--'); hold on

drawArrow([0 x0],[0 y0],xran,yran,'linewidth',2); hold on

%%% rotate the coordinate system by 45 deg, clockwise, with rotation
%%% matrix, see the expression in new coords
ang = deg2rad(45);
R = [cos(ang) -sin(ang); sin(ang) cos(ang)];
newcoord = R * [x0; y0];    % matrix multiplication
x1 = newcoord(1)
y1 = newcoord(2)
drawArrow([0 4],[0 4],xran,yran,'color',[128/255 128/255 128/255],'linewidth',1); hold on
% drawArrow([0 x1],[0 y1],xran,yran,'color','b','linewidth',1.5); hold on


%%% rotate the record by 15 deg, clockwise, with complex rotation.
%%% expression in original system.
c = x0+1i*y0;
ang = deg2rad(15);
rotc = c * exp(-1i * ang);
x2 = real(rotc)
y2 = imag(rotc)
drawArrow([0 x2],[0 y2],xran,yran,'color','r','linewidth',1.5); hold on


%%% rotate the coordinate system by 15 deg, counter-clockwise, with
%%% rotation matrix,
%%% should be exactly the same as  
ang = deg2rad(-15);
R = [cos(ang) -sin(ang); sin(ang) cos(ang)];
newcoord = R * [x0; y0];    % matrix multiplication
x3 = newcoord(1)
y3 = newcoord(2)
drawArrow([0 x3],[0 y3],xran,yran,'color','g','linewidth',1.5); hold on

