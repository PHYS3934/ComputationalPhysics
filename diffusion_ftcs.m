function diffusion_ftcs(tau,h)
% Solve the 1-D diffusion equation for an initial spike profile
% with Dirichlet conditions using FTCS, in a matrix formulation

% Clear memory and show only a few digits
% clear('all');
% format('short');

% Thermal conductivity
kappa = 1;

% Time step and spatial step
if nargin < 1
    tau = 1e-4;
end
if nargin < 2
    h = 0.02;
end

% Number of time steps
numSteps = 50;
frameUpdateLag = 0.2;
zoomIn = false;

% Calculate the ratio tau/(th/2), where th is the approximate diffusion
% time for one spatial step h
th = h^2/kappa;
disp(['Ratio tau/(0.5*th): ',num2str(tau/(0.5*th))]);

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

% Initial conditions, temp0: a spike at x = 1/2
temp0 = zeros(L,1);
temp0(round(L/2)) = 1/h;
temp = temp0; % temp updates across the run

% Record T(x,t) matrix for visualization
time = tau*(0:numSteps);
temp_xt = zeros(L,numSteps+1);
temp_xt(:,1) = temp;

%-------------------------------------------------------------------------------
% Plot initial condition and set up for animation
f = figure('color','w');
hold('on')
niceRed = [0.84,0.09,0.11];
niceOrange = [0.99,0.68,0.38];
niceBlue = [0.17,0.51,0.73];
temp_an = zeros(L,1); % (dummy)
p_TempAnal = plot(x,temp_an,'k','LineWidth',1.5); % analytical profile
p_Temp0 = plot(x,temp0,'-','Color',niceBlue,'LineWidth',1.5); % initial profile
p_Temp = plot(x,temp,'o-','Color',niceRed,...
                'MarkerFaceColor',niceRed,...
                'MarkerEdgeColor',niceOrange);
h_legend = legend('Analytic solution','Initial profile','Numerical solution');
h_legend.Box = 'off';
h_legend.Location = 'NorthWest';
xlabel('Position (non-dim.)');
ylabel('Temperature (non-dim.)');

%-------------------------------------------------------------------------------
% March forwards in time, FTCS style!
for n = 1:numSteps

    % Update the temperature profile
    temp = A*temp;

    % Recalculate the profile for the (approximate) analytic solution
    sig = sqrt(2*kappa*time(n+1));
    temp_an = exp(-(x - 0.5).^2/(2*sig^2))/(sqrt(2*pi)*sig);

    % Animation:
    title(sprintf('Time: %.2g (%u/%u)',time(n+1),n,numSteps));
    p_Temp.YData = temp; % update current profile
    p_TempAnal.YData = temp_an; % update analytical profile
    % Follow evolution more closely near end:
    if zoomIn && (n > numSteps/2)
        axis([0 1 min(temp_an) max(temp_an)]);
    end
    drawnow()
    pause(frameUpdateLag);

    % Record T(x,t) matrix for visualization
    temp_xt(:,n+1) = temp;
end

%-------------------------------------------------------------------------------
% Visualization of temperature versus position and time.
% Matlab uses a surfc(x(i),y(j),z(j,i)) convention, hence the transpose.
showSurfacePlot = false;
if showSurfacePlot
    figure(2);
    surfc(x,time,temp_xt');
    shading('flat');
    xlabel('Position'); ylabel('Time'); zlabel('Temperature');
    colormap('hot')
end
