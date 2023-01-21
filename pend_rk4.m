% pend_rk4.m
% Solve the nonlinear pendulum problem using RK4

% Clear memory and show only a few digits
clear('all');
format('short');

%-------------------------------------------------------------------------------
% Time grid parameters:
tau = 0.05; % Time step
T = 8; % Total integration time
numSteps = ceil(T/tau); % Number of time steps
frameUpdateLag = 0; % Slow down the animation

%-------------------------------------------------------------------------------
% Initial conditions
theta1 = 50; % Initial angle in degrees
theta = theta1*pi/180;
x = [theta 0];

%-------------------------------------------------------------------------------
% Rainbow Mode
rainbowMode = true;
rainbow = jet(64);

%-------------------------------------------------------------------------------
% Preallocate:
time = (0:numSteps)*tau;

%-------------------------------------------------------------------------------
% Set up plot for animation:
figure('color','w');
h = plot([0,sin(x(1))],[0,-cos(x(1))],'o-k','LineWidth',2);
axis('equal')
axis([-1 1 -1 1])
% Prevent axes limits from changing:
h.YLimInclude = 'off';
h.XLimInclude = 'off';

%-------------------------------------------------------------------------------
% Fourth-order Runge-Kutta integration
for n = 1:numSteps

    % One step of RK4
    f1 = rhs_pend(x);
    f2 = rhs_pend(x + 0.5*tau*f1);
    f3 = rhs_pend(x + 0.5*tau*f2);
    f4 = rhs_pend(x + tau*f3);
    x = x + tau*(f1 + 2*f2+2*f3 + f4)/6;

    % Update the pendulum position:
    % Co-ordinates of the pendulum bar
    xPend = [0  sin(x(1))];
    yPend = [0 -cos(x(1))];
    % plot(xPend,yPend,'o-'):
    h.XData = xPend;
    h.YData = yPend;
    if rainbowMode
        h.Color = rainbow(rem(n,64) + 1,:);
    end

    % Update time and keep axes fixed:
    title(sprintf('Time: %5.3f',time(n) + tau))
    pause(frameUpdateLag);
    drawnow()
end
