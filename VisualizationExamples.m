%% Create Sample data for visualization

%% Electrode 1 (Context 1 + Context 2)

m=1;
N=60*4;
sig1=1;
sig2=0.5;

mean1=0;
mean2=2.5;

f1=1/45;
f2=1/3;
Amp1=2.5;
Amp2=1.25;
ts=1/4;
T=60;
t=0:ts:T-1/T;
slow=Amp1*sin(2*pi*f1*t);
fast=Amp2*sin(2*pi*f2*t);
%y = pinknoise(m, N);
e1c1=randn(m,N)*sqrt(sig1)+mean1+slow+0.25*fast;
%y = pinknoise(m, N);
e1c2=randn(m,N)*sqrt(sig2)+mean2+fast+0.15*slow;

%y = pinknoise(m, N);
e2c1=randn(m,N)*sqrt(sig1)+mean1+slow+0.25*fast;
%y = pinknoise(m, N);
e2c2=randn(m,N)*sqrt(sig2)+mean2+fast*-1+0.15*slow;


subplot(1,2,1)
plot(t,e1c1, 'b', 'LineWidth', 1.5);
hold on;
plot(t,e2c1, 'Color', [255 69 0]./255, 'LineWidth', 1.5);
xlim([0 60]);
ylim([-6 6]);
title('State  1');
legend('Electrode 1', 'Electrode 2');
% Electrode 1 (Context 1 + Context 2)
xlabel('Time (seconds)');
ylabel('Amplitude a.u');


subplot(1,2,2)
plot(t,e1c2, 'b', 'LineWidth', 1.5);
hold on;
plot(t,e2c2, 'Color', [255 69 0]./255, 'LineWidth', 1.5);
xlim([0 60]);
title('State  2');
legend('Electrode 1', 'Electrode 2');

ylim([-6 6]);
xlabel('Time (seconds)');
ylabel('Amplitude a.u');
%%

figure
%subplot(1,3,3)
scatter(e1c1, e2c1, 30, 'MarkerFaceColor', [0.466 0.6740 0.1880], 'MarkerEdgeColor', [0.466 0.6740 0.1880]);
xlabel('Electrode 1');
ylabel('Electrode 2');
%title('State 1');
%subplot(2,2,4)

hold on;
scatter(e1c2, e2c2,  30,'MarkerFaceColor', [0.494 0.184 0.556], 'MarkerEdgeColor',  [0.466 0.6740 0.1880]);

%%


data=[e1c1; e2c1]';
% Calculate the eigenvectors and eigenvalues
covariance = cov(data);
[eigenvec, eigenval ] = eig(covariance);

% Get the index of the largest eigenvector
[largest_eigenvec_ind_c, r] = find(eigenval == max(max(eigenval)));
largest_eigenvec = eigenvec(:, largest_eigenvec_ind_c);

% Get the largest eigenvalue
largest_eigenval = max(max(eigenval));

% Get the smallest eigenvector and eigenvalue
if(largest_eigenvec_ind_c == 1)
    smallest_eigenval = max(eigenval(:,2))
    smallest_eigenvec = eigenvec(:,2);
else
    smallest_eigenval = max(eigenval(:,1))
    smallest_eigenvec = eigenvec(1,:);
end

% Calculate the angle between the x-axis and the largest eigenvector
angle = atan2(largest_eigenvec(2), largest_eigenvec(1));

% This angle is between -pi and pi.
% Let's shift it such that the angle is between 0 and 2pi
if(angle < 0)
    angle = angle + 2*pi;
end

% Get the coordinates of the data mean
avg = mean(data);

% Get the 95% confidence interval error ellipse
chisquare_val = 2.4477;
theta_grid = linspace(0,2*pi);
phi = angle;
X0=avg(1);
Y0=avg(2);
a=chisquare_val*sqrt(largest_eigenval);
b=chisquare_val*sqrt(smallest_eigenval);

% the ellipse in x and y coordinates 
ellipse_x_r  = a*cos( theta_grid );
ellipse_y_r  = b*sin( theta_grid );

%Define a rotation matrix
R = [ cos(phi) sin(phi); -sin(phi) cos(phi) ];

%let's rotate the ellipse to some angle phi
r_ellipse = [ellipse_x_r;ellipse_y_r]' * R;

% Draw the error ellipse
plot(r_ellipse(:,1) + X0,r_ellipse(:,2) + Y0,'-', 'color', 'g')
hold on;


% Plot the eigenvectors
%quiver(X0, Y0, largest_eigenvec(1)*sqrt(largest_eigenval), largest_eigenvec(2)*sqrt(largest_eigenval), '-m', 'LineWidth',2);
%quiver(X0, Y0, smallest_eigenvec(1)*sqrt(smallest_eigenval), smallest_eigenvec(2)*sqrt(smallest_eigenval), '-g', 'LineWidth',2);
%hold on;

data=[e1c2; e2c2]';
% Calculate the eigenvectors and eigenvalues
covariance = cov(data);
[eigenvec, eigenval ] = eig(covariance);

% Get the index of the largest eigenvector
[largest_eigenvec_ind_c, r] = find(eigenval == max(max(eigenval)));
largest_eigenvec = eigenvec(:, largest_eigenvec_ind_c);

% Get the largest eigenvalue
largest_eigenval = max(max(eigenval));

% Get the smallest eigenvector and eigenvalue
if(largest_eigenvec_ind_c == 1)
    smallest_eigenval = max(eigenval(:,2))
    smallest_eigenvec = eigenvec(:,2);
else
    smallest_eigenval = max(eigenval(:,1))
    smallest_eigenvec = eigenvec(1,:);
end

% Calculate the angle between the x-axis and the largest eigenvector
angle = atan2(largest_eigenvec(2), largest_eigenvec(1));

% This angle is between -pi and pi.
% Let's shift it such that the angle is between 0 and 2pi
if(angle < 0)
    angle = angle + 2*pi;
end

% Get the coordinates of the data mean
avg = mean(data);

% Get the 95% confidence interval error ellipse
chisquare_val = 2.4477;
theta_grid = linspace(0,2*pi);
phi = angle;
X0=avg(1);
Y0=avg(2);
a=chisquare_val*sqrt(largest_eigenval);
b=chisquare_val*sqrt(smallest_eigenval);

% the ellipse in x and y coordinates 
ellipse_x_r  = a*cos( theta_grid );
ellipse_y_r  = b*sin( theta_grid );

%Define a rotation matrix
R = [ cos(phi) sin(phi); -sin(phi) cos(phi) ];

%let's rotate the ellipse to some angle phi
r_ellipse = [ellipse_x_r;ellipse_y_r]' * R;

% Draw the error ellipse
plot(r_ellipse(:,1) + X0,r_ellipse(:,2) + Y0,'-', 'color', 'r')
hold on;

% Plot the eigenvectors
%quiver(X0, Y0, largest_eigenvec(1)*sqrt(largest_eigenval), largest_eigenvec(2)*sqrt(largest_eigenval), '-m', 'LineWidth',2);
%quiver(X0, Y0, smallest_eigenvec(1)*sqrt(smallest_eigenval), smallest_eigenvec(2)*sqrt(smallest_eigenval), '-g', 'LineWidth',2);
%hold on;




% Set the axis labels
xlabel('Electrode 1');
ylabel('Electrode 2');

legend('State 1', 'State 2')

%title('State 2');