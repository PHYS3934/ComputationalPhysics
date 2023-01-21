function diffusion_ftcs_leftRightBoundary(f,leftOrRight)
% Solve the 1-D diffusion equation for some initial profile
% with Dirichlet conditions using FTCS, using a matrix formulation.

%-------------------------------------------------------------------------------
% INPUTS:
%-------------------------------------------------------------------------------
% Set the stability factor
if nargin < 1
    f = 0.25;
end

% leftOrRight
if nargin < 2
    % leftOrRight = 'both'; % Use this for Q6
    leftOrRight = 'right'; % Use this for Q8
end
%-------------------------------------------------------------------------------

% Turn on to plot snapshots
plotSnapshots = true;

% Thermal conductivity
kappa = 5;

% Spatial step
h = 0.05;

% Time step
tau = f*h^2/kappa; % Time step

% Number of time steps
numSteps = 400;
frameUpdateLag = 0;

% Column vector of x values
x = (0:h:1)';
L = length(x);

% Construct the matrix D associated with the second spatial
% derivative and the boundary conditions
D = -2*eye(L);
D = D + diag(ones(L-1,1),+1) + diag(ones(L-1,1),-1);
D = kappa*tau*D/h^2;

% Impose the Dirichlet boundary conditions
D(1,:) = zeros(1,L);
D(L,:) = zeros(1,L);

% Construct the update matrix
A = eye(L) + D;

% Initial conditions: increased temperature at ends (indices 1 and L).
temp0 = zeros(L,1);
switch leftOrRight
case 'both'
    temp0(1) = 1; % KEEP FOR Q6, COMMENT OUT FOR Q8.
    temp0(L) = 1;
case 'right'
    temp0(L) = 1;
case 'left'
    temp0(1) = 1;
end
temp = temp0; % temp updates across the run.

% Record T(x,t) matrix for visualisation
time = tau*(0:numSteps);

%-------------------------------------------------------------------------------
% Plot initial condition and set up for animation
f = figure('color','w');
hold('on')
temp_an = zeros(L,1); % (dummy)
niceBlue = [0.17,0.51,0.73];
niceRed = [0.84,0.09,0.11];
niceOrange = [0.99,0.68,0.38];
p_Temp0 = plot(x,temp0,'-','Color',niceBlue,'LineWidth',1.5); % initial profile
p_Temp = plot(x,temp,'o-','Color',niceRed,...
                'MarkerFaceColor',niceRed,...
                'MarkerEdgeColor',niceOrange);
xlabel('Position (non-dim.)');
ylabel('Temperature (non-dim.)');

%-------------------------------------------------------------------------------
% March forwards in time, FTCS style!
for n = 1:numSteps

    % Update the temperature profile
    temp = A*temp;

    % Animation:
    title(sprintf('Time: %.2g (%u/%u)',time(n+1),n,numSteps));
    p_Temp.YData = temp; % update current profile

    % Plot regular snapshots:
    if plotSnapshots && rem(n,50)==0
        plot(x,temp,'o-','Color',niceRed,...
                'MarkerFaceColor',niceRed,...
                'MarkerEdgeColor',niceOrange);
    end

    drawnow()
    pause(frameUpdateLag);
end

legend([p_Temp0,p_Temp],{'Initial Profile','Final temperature profile'},'Location','northwest')
