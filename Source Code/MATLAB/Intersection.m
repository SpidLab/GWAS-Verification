xaxis = dlmread('index_sample.txt'); % read the x-axis values
correct_y = dlmread('correct_sample.txt'); % read the correct values which are used as y coordinate
incorrect_y = dlmread('incorrect_sample.txt');  % read the correct values which are used as y coordinate

% xaxis = dlmread('index.txt'); % read the x-axis values
% correct_y = dlmread('correct.txt'); % read the correct values which are used as y coordinate
% incorrect_y = dlmread('incorrect.txt');  % read the correct values which are used as y coordinate

x1=xaxis(1,:)
y1=correct_y(1,:)
x2=xaxis(1,:);
y2=incorrect_y(1,:);
P = InterX([x1;y1],[x2;y2]);%compute the intersection points
plot(x1,y1,x2,y2,P(1,:),P(2,:),'ro') %plot the correct and incorrect lines
disp(P) %print the intersection points