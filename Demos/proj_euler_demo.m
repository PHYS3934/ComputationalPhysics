function [range_m,an_range_m] = proj_euler(speed_m,angle,tau)
%-------------------------------------------------------------------------------
% Simple projectile motion using Euler's method
%-------------------------------------------------------------------------------

% Dimensionalisation parameters
G = 9.8; % Acceleration due to gravity (m/s^2)
Ls = 1.0; % Choice for scaling length (m)
Ts = sqrt(Ls/G); % Choice for scale for time (s)

% Prompt user for initial speed and angle
% speed_m = input('Enter initial speed in m/s: ');
% angle = input('Enter initial angle in degrees: ');

% Convert angle to radians
angle = angle*pi/180;

% Non-dimensionalise initial speed
speed = speed_m/(Ls/Ts);

% Row vectors for non-dimensional position and velocity
pos = [0 0];
vel = speed*[cos(angle) sin(angle)];

% Store Initial Condition (for plotting):
x(1) = pos(1);
y(1) = pos(2);

%-------------------------------------------------------------------------------
% Euler's method!:
%-------------------------------------------------------------------------------
f = figure('color','w');
hold('on')
xlabel('Distance (m)')
ylabel('Height (m)')
n = 1; % Index
while pos(2) >= 0
    % Compute one step of Euler's method:
    % r_{n+1} = r_n + τ v_n
    pos = pos + tau*vel;
    % v_{n+1} = v_n - τ \hat{y}
    vel = vel + tau*[0 -1];

    % Update the index:
    n = n + 1;

    % Store position for plotting:
    x(n) = pos(1);
    y(n) = pos(2);

    plot(Ls*x(n-1:n),Ls*y(n-1:n),'o-k')
    drawnow();
    pause(0.01)
end

%-------------------------------------------------------------------------------
% Estimate range
%-------------------------------------------------------------------------------
% Linear interpolation to estimate the range of the projectile
coOrdsOver = [x(end-1),y(end-1)]; % last point projectile above axis
coOrdsUnder = [x(end),y(end)]; % projectile under ground

range = coOrdsUnder(1) - coOrdsUnder(2)*(coOrdsUnder(1)-coOrdsOver(1))/(coOrdsUnder(2)-coOrdsOver(2));
range_m = Ls*range; % convert back to m

% Analytic expression for range
an_range_m = speed_m^2*sin(2*angle)/G;

%-------------------------------------------------------------------------------
% Plot the trajectory (as dimensional values):
%-------------------------------------------------------------------------------
plot(range_m,0,'ok','MarkerFaceColor','r')
plot(an_range_m,0,'ok','MarkerFaceColor','g')
% Force minimum vertical axis limit to zero:
ax = gca();
ax.YLim(1) = 0; %


% Put information out to command line:
fprintf('Range (m): %f\n',range_m)
fprintf('Analytic value for range (m): %f.\n',an_range_m)

end
