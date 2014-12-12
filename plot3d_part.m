function plot3d_part(thislag,varargin)

spacing=1;
title_str='';


if nargin>=1
    spacing=varargin{1};
end
if nargin>=2
    title_str=varargin{2};
end



[a b]=size(thislag.x);

figure

subplot(2,2,1)
hold on
for i=1:spacing:a
plot3(thislag.x(i,:),thislag.y(i,:),thislag.z(i,:),'r')
plot3(thislag.x(i,:),thislag.y(i,:),thislag.z(i,:),'.g')
plot3(thislag.x(i,1),thislag.y(i,1),thislag.z(i,1),'.k')
end
xlabel('x')
ylabel('y')
zlabel('z')
view(45,25)



subplot(2,2,2)

hold on
for i=1:spacing:a
plot(thislag.x(i,:),thislag.y(i,:),'r')
plot(thislag.x(i,:),thislag.y(i,:),'.g')
plot(thislag.x(i,1),thislag.y(i,1),'.k')
xlabel('x')
ylabel('y')
end



subplot(2,2,3)

hold on
for i=1:spacing:a
plot(thislag.x(i,:),thislag.z(i,:),'r')
plot(thislag.x(i,:),thislag.z(i,:),'.g')
plot(thislag.x(i,1),thislag.z(i,1),'.k')
xlabel('x')
ylabel('z')
end



subplot(2,2,4)

hold on
for i=1:spacing:a
plot(thislag.y(i,:),thislag.z(i,:),'r')
plot(thislag.y(i,:),thislag.z(i,:),'.g')
plot(thislag.y(i,1),thislag.z(i,1),'.k')
xlabel('y')
ylabel('z')
end


suptitle(title_str)
