% kepler_euler.m
% Motion under a central force, using Euler's method
%-------------------------------------------------------------------------------

% Clear memory and only show a few digits
clear('all');
format('short');

% Time step (non-dim.)
tau = 0.05;

% Initial position (non-dim.) - this should be fixed
pos = [1 0];

% Initial velocity (non-dim.) - vary the y-component
vel = [0 1];

% Total integration time
T = 4*pi;

% Number of integration steps
numSteps = ceil(T/tau);

% Plot only 100 frames in total
numFrames = 100;
% Plot every 'skip' iterations:
skip = ceil(numSteps/numFrames);

% Calculate trajectory from analytic solution.
% See Appendix of Lecture 2.
[xan,yan] = kepler_analytic(vel,T);

% Preallocate vectors for speed:
x = zeros(numSteps,1);
y = zeros(numSteps,1);
time = zeros(numSteps,1);
energy = zeros(numSteps,1);

%-------------------------------------------------------------------------------
% Euler's method integration
%-------------------------------------------------------------------------------
figure(1);
xlabel('x'); ylabel('y');
for n = 1:numSteps

    % Store position and time for plotting: time step n
    x(n) = pos(1);
    y(n) = pos(2);
    time(n) = (n-1)*tau;

    % Plot numerical and analytic solution
    if rem(n,skip)==0
        plot(x,y,'g-',pos(1),pos(2),'ko',xan,yan,'b',0,0,'ro')
        title(sprintf('Time: %f',time(n)));
        axis('equal'); % Preserve aspect ratio
        drawnow; % Draw immediately
    end

    % Calculate radial position, speed and acceleration at step n
    r = norm(pos);
    speed = norm(vel);
    accel = -pos/r^3;

    % Calculate total energy at step n and store
    energy(n) = 0.5*speed^2-1/r;

    % One step of Euler's method
    pos = pos + tau*vel;
    vel = vel + tau*accel;
end

% Plot energy versus time
figure(2);
plot(time,energy);
xlabel('Time (non-dim.)');
ylabel('Total energy (non-dim.)');
