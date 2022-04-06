% proj_euler
% ---------------------------------------------
% Simple projectile motion using Euler's method

% Only show a few digits:
format('short');

% Dimensionalisation parameters
G = 9.8; % Acceleration due to gravity (m/s^2)
Ls = 1.0; % Choice for scaling length (m)
Ts = sqrt(Ls/G); % Choice for scale for time (s)

% Non-dimensional time-step, tau
tau = 0.5;

% Prompt user for initial speed and angle
speed_m = input('Enter initial speed in m/s: ');
angle = input('Enter initial angle in degrees: ');

% Convert angle to radians
angle = angle*pi/180;

% Non-dimensionalise initial speed
speed = speed_m/(Ls/Ts);

% Row vectors for non-dimensional position and velocity
pos = [0 0];
vel = speed*[cos(angle) sin(angle)];

% Initialize variables to save for plotting:
x = [];
y = [];

%-------------------------------------------------------------------------------
% Euler's method!:
%-------------------------------------------------------------------------------
n = 0;
while pos(2) >= 0
    n = n + 1;

    % Store position for plotting:
    x(n) = pos(1);
    y(n) = pos(2);

    % Compute one step of Euler's method:
    % r_{n+1} = r_n + τ v_n
    pos = pos + tau*vel;
    % v_{n+1} = v_n - τ \hat{y}
    vel = vel + tau*[0 -1];
end
%-------------------------------------------------------------------------------

% Plot the trajectory (as dimensional values):
f = figure('color','w');
plot(Ls*x,Ls*y,'o-k')
xlabel('Distance (m)')
ylabel('Height (m)')

% Linear interpolation to estimate the range of the projectile
range = pos(1) - pos(2)*(pos(1)-x(end))/(pos(2)-y(end));
range_m = Ls*range; % convert back to m
fprintf('Range (m): %f\n',range_m)

% Analytic expression for range:
an_range_m = speed_m^2*sin(2*angle)/G;
fprintf('Analytic value for range (m): %f.\n',an_range_m)
