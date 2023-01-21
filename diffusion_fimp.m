function diffusion_fimp(fac)
% Solve the 1-D diffusion equation using the fully implicit scheme
% with Dirichlet BCs

% Defaults:
if nargin < 1
    fac = 20; % For FTCS, this must be <= 1/2 for stability
end

% Fixed parameters
kappa = 1;           % Thermal conductivity
h = 0.01;            % Spatial step
numSteps = 200;      % Number of time steps
tau = fac*h^2/kappa; % Time step
snapshotPeriod = 20; % Take regular snapshots (set to zero for none)
showHeatmap = false; % Show the heatmap plot at the end.

% Display the FTCS stability factor
fprintf('FTCS stability factor: %g (FTCS stable if <1/2)\n',fac);

%-------------------------------------------------------------------------------
% Column vector of x values
x = (0:h:1)';
L = length(x);

% Construct the matrix D associated with the second spatial derivative
D = -2*eye(L);
D = D + diag(ones(L-1,1),+1) + diag(ones(L-1,1),-1);

% Construct the matrix in the linear system being solved at each step
M = eye(L) - fac*D;

% Impose Dirichlet boundary conditions
M(1,:) = zeros(1,L);
M(1,1) = 1;
M(L,:) = zeros(1,L);
M(L,L) = 1;

% Initial conditions and BCs
temp = zeros(L,1);
temp(1) = 0;
temp(L) = 1;

% Record T(x,t) matrix for visualisation
temp_xt = zeros(L,numSteps+1);
temp_xt(:,1) = temp;
time = tau*(0:numSteps);

%-------------------------------------------------------------------------------
% Plot initial condition and set up for animation
f = figure('color','w');
hold('on')
niceRed = [0.84,0.09,0.11];
niceOrange = [0.99,0.68,0.38];
niceBlue = [0.17,0.51,0.73];
p_TempEqm = plot(x,x,'-','Color',niceBlue,'LineWidth',1.5); % initial profile
p_Temp = plot(x,temp,'o-','Color',niceRed,...
                'MarkerFaceColor',niceRed,...
                'MarkerEdgeColor',niceOrange);
h_legend = legend([p_TempEqm,p_Temp],{'Equilibrium','Numerical solution'});
h_legend.Box = 'off';
h_legend.Location = 'NorthWest';
xlabel('Position (non-dim.)');
ylabel('Temperature (non-dim.)');
axis([0 1 0 1]);

%-------------------------------------------------------------------------------
% March forwards in time
for n = 1:numSteps

    % Solve the matrix equation for a time step
    temp = M\temp;

    % Plot the temperature versus position and the equilibrium
    % solution, with annotations
    p_Temp.YData = temp;
    title(sprintf('f = %g, Time: %g [%u/%u]',fac,time(n+1),n,numSteps));
    drawnow()

    if snapshotPeriod > 0 && rem(n,snapshotPeriod)==0
        plot(x,temp,'o-','Color',niceRed,...
                'MarkerEdgeColor',niceRed,...
                'MarkerFaceColor',niceRed,...
                'MarkerSize',3);
    end

    % Record T(x,t) matrix for visualisation
    temp_xt(:,n+1) = temp;
end

%-------------------------------------------------------------------------------
% Visualize temperature versus position and time as a heatmap
if showHeatmap
    f2 = figure(2);
    f2.Color = 'w';
    colormap(hot)
    imagesc(x,time,flipud(temp_xt'));
    xlabel('Position');
    ylabel('Time');
    cB = colorbar();
    cB.Label.String = 'Temperature';
end

end
