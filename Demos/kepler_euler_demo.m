function kepler_euler_demo(tau,T)
% Motion under a central force, using Euler's method
%-------------------------------------------------------------------------------
% Time step (non-dim.)
if nargin < 1
    tau = 0.05;
end

% Total integration time
if nargin < 2
    T = 4*pi;
end

% Initial position (non-dim.) - this should be fixed
pos = [1 0];

% Initial velocity (non-dim.) - vary the y-component
vel = [0 1];

% Number of integration steps
numSteps = ceil(T/tau);

% Plot only a subset of frames
numFrames = 50; % numSteps, 100
% Plot every 'skip' iterations:
skip = ceil(numSteps/numFrames);

% Calculate trajectory from analytic solution.
% See Appendix of Lecture 2.
[xan,yan] = kepler_analytic(vel,T);

% Preallocate vectors for speed:
time = tau*(0:numSteps);
x = zeros(numSteps+1,1);
y = zeros(numSteps+1,1);
energy = zeros(numSteps+1,1);

% Initial values:
x(1) = pos(1);
y(1) = pos(2);
r = norm(pos);
speed = norm(vel);
accel = -pos/r^3;
energy(1) = 0.5*speed^2 - 1/r;

%-------------------------------------------------------------------------------
% Set up figure:
niceRed = [0.84,0.09,0.11];
niceOrange = [0.99,0.68,0.38];
niceBlue = [0.17,0.51,0.73];
figure('color','w');
subplot(1,2,1)
plot(xan,yan,'color',niceBlue)
hold('on')
plot(0,0,'o','color',niceBlue)
xlabel('x'); ylabel('y');
subplot(1,2,2)
hold('on')
xlabel('Time (non-dim.)');
ylabel('Total energy (non-dim.)');

%-------------------------------------------------------------------------------
% Euler's method integration
%-------------------------------------------------------------------------------
for n = 1:numSteps

    % Plot numerical and analytic solution
    if n > 1 && rem(n,skip)==0
        subplot(1,2,1)
        plot(x(n-1:n),y(n-1:n),'-','color',niceRed);
        plot(x(n),y(n),'.','color',niceRed)
        title(sprintf('Time: %f',time(n)));
        axis('equal'); % Preserve aspect ratio
        subplot(1,2,2)
        plot(time(n),energy(n),'.k');
        drawnow(); % Draw immediately
        % pause(0.01)
    end

    % Take one step of Euler's method:
    pos = pos + tau*vel;
    vel = vel + tau*accel;

    % Calculate radial position, speed and acceleration at step n:
    r = norm(pos);
    speed = norm(vel);
    accel = -pos/r^3;

    % Store position, time, energy after taking the step:
    x(n+1) = pos(1);
    y(n+1) = pos(2);
    energy(n+1) = 0.5*speed^2 - 1/r;
end
