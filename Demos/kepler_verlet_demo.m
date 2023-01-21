function kepler_verlet_demo(tau,T)
% Motion under a central force, using the Verlet method
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

% Plot only a subset of frames in total
numFrames = 100;
% Plot every 'skip' iterations:
skip = ceil(numSteps/numFrames);

% Calculate trajectory from analytic solution (see Appendix of Lecture 2).
[xan,yan] = kepler_analytic(vel,T);

% Preallocate vectors for speed:
time = tau*(0:numSteps);
x = zeros(numSteps+1,1);
y = zeros(numSteps+1,1);
energy = zeros(numSteps+1,1);

% Set initial values:
x(1) = pos(1);
y(1) = pos(2);
speed = norm(vel);
r = norm(pos);
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
% Verlet method integration
%-------------------------------------------------------------------------------
for n = 1:numSteps

    % Plot numerical and analytic solution
    if rem(n,skip)==0 && n > 1
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
    %
    %     plot(x(1:n),y(1:n),'g-',x(n),y(n),'ko',xan,yan,'b',0,0,'ro')
    %     title(sprintf('Time: %f',time(n)));
    %     axis('equal'); % Preserve aspect ratio
    %     drawnow(); % Draw immediately
    % end

    % Take one step of the Verlet Method to update position
    if n==1
        % Get started with a midpoint method step
        next = pos + tau*vel + 0.5*tau^2*accel;
    else
        % Normal Verlet update:
        next = 2*pos - prev + tau^2*accel;
    end

    % Take one step of the Verlet Method to update velocity
    % (only required to compute energy, not the trajectory)
    if n==1
        % Do a special (Velocity-Verlet) update for the very first step:
        r = norm(next);
        accel_n2 = -next/r^3;
        vel = vel + tau/2*(accel + accel_n2);
    else % n > 1
        % Do a Verlet update:
        vel = (next - prev)/(2*tau);
    end

    % Calculate speed, radial position, and acceleration after taking the step:
    speed = norm(vel);
    r = norm(next);
    accel = -next/r^3;

    % Update energy after taking the step:
    energy(n+1) = 0.5*speed^2 - 1/r;

    % Store position update:
    x(n+1) = next(1);
    y(n+1) = next(2);

    % Update 'prev' and 'pos' to calculate 'next' in the following step:
    prev = pos;
    pos = next;
end

%-------------------------------------------------------------------------------
% Plot energy versus time
% figure(2);
% hold('on')
% plot(time,energy);
% plot([min(time),max(time)],ones(2,1)*energy(1),'--k')
% xlabel('Time (non-dim.)');
% ylabel('Total energy (non-dim.)');
