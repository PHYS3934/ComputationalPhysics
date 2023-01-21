function diffusion_cnic(tau,h)
% Solve the diffusion equation using Crank-Nicolson
%-------------------------------------------------------------------------------

if nargin < 1
    % time step
    tau = 8e-4;
end
if nargin < 2
    % spatial step
    h = 0.02;
end

% Thermal conductivity
kappa = 1;

% Number of time steps
totalTime = 0.02;
numSteps = ceil(totalTime/tau);
frameUpdateLag = 0.1; % pause between updates
zoomIn = true;
showHeatMap = false;

% Display value of FTCS stability factor
fac = kappa*tau/h^2;
disp(['FTCS stability factor (stable <= 1/2): ',num2str(fac)]);

% Vector of x values
x = (0:h:1)';
L = length(x);

% Construct the matrix D associated with the second spatial derivative
D = -2*eye(L);
D = D + diag(ones(L-1,1),+1) + diag(ones(L-1,1),-1);

% Impose Dirichlet boundary conditions
D(1,:) = zeros(1,L);
D(L,:) = zeros(1,L);

% Construct the matrix for the linear system solved at each step of
% Crank-Nicolson
A = 0.5*(eye(L) - 0.5*fac*D);

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
figure('color','w');
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
h_legend = legend('Approx analytic solution','Initial profile','Numerical solution');
h_legend.Box = 'off';
h_legend.Location = 'NorthWest';
xlabel('Position (non-dim.)');
ylabel('Temperature (non-dim.)');

%-------------------------------------------------------------------------------
% March forwards in time, Crank-Nicholson Style!
for n = 1:numSteps

    % Perform a Crank-Nicolson update
    chi = A\temp;
    temp = chi - temp;

    % Recalculate the profile for the (approximate) analytic solution
    sig = sqrt(2*kappa*time(n+1));
    temp_an = exp(-(x - 0.5).^2/(2*sig^2))/(sqrt(2*pi)*sig);

    % Animation:
    title(sprintf('f = %g; Time: %.2g (%u/%u)',...
                        kappa*tau/h^2,time(n+1),n,numSteps));
    p_Temp.YData = temp; % update current profile
    p_TempAnal.YData = temp_an; % update analytical profile
    if zoomIn
        ylim([0,max(temp)*1.5])
    end
    drawnow()
    pause(frameUpdateLag);

    % Record T(x,t) matrix for visualization
    temp_xt(:,n+1) = temp;
end

%-------------------------------------------------------------------------------
%-------------------------------------------------------------------------------
% Visualize temperature versus position and time. Matlab
% uses a surfc(x(i),y(j),z(j,i)) convention, hence the transpose
if showHeatMap
    f2 = figure(2);
    f2.Color = 'w';
    ax = gca();
    colormap(hot)
    imagesc(x,time,temp_xt');
    ax.YDir = 'normal';
    xlabel('Position');
    ylabel('Time');
    cB = colorbar();
    cB.Label.String = 'Temperature';
end

end
