function createfigure(X1, Y1)
%CREATEFIGURE(X1, Y1)
%  X1:  vector of x data
%  Y1:  vector of y data

%  Auto-generated by MATLAB on 18-Feb-2019 04:50:40

% Create figure
figure('OuterPosition',[672 538 576 513]);

% Create axes
axes1 = axes;
hold(axes1,'on');

% Create plot
plot(X1,Y1);

% Create ylabel
ylabel('Weight [lbs]');

% Create xlabel
xlabel('Range [mi]');

% Create title
title('Range Trade Study');

box(axes1,'on');
