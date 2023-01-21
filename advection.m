function advection(whatMethod,tau)
% advection Solve the advection equation with periodic BCs

%% Set defaults:
% Set default method:
if nargin < 1
    whatMethod = 'ftcs';
    % whatMethod = 'lax';
end
fprintf(1,'Using the %s method!\n',whatMethod);
% Set default time step, tau:
if nargin < 2
    tau = 0.01;
end

%--------------------------------------------------------------------------
%% Set parameters:
%--------------------------------------------------------------------------
% Propagation speed:
c = 1;

% Spatial step:
h = 0.01;

% Lag time between plot updates
frameUpdateLag = 0.05; % (s)

%--------------------------------------------------------------------------
%% Initialization
%--------------------------------------------------------------------------

% Number of steps such that the pulse should propagate through
% the periodic domain once
numSteps = ceil((1 + h)/(c*tau));

% Vector of x values
x = 0:h:1;
L = length(x);

%-------------------------------------------------------------------------------
%% Construct the update matrix, M
%-------------------------------------------------------------------------------
g = c*tau/h;
switch whatMethod
case 'ftcs'
    M = construct_update_matrix(L,'ftcs_advection',g);
case 'lax'
    M = construct_update_matrix(L,'lax_advection',g);
otherwise
    error('Unknown method, ''%s''',whatMethod)
end

% Calculate the spectral radius of M
rho = max(abs(eig(M)));
disp(['Spectral radius: ',num2str(rho)]);

%-------------------------------------------------------------------------------
%% Set Initial conditions (a Gaussian pulse at x = 0.5)
%-------------------------------------------------------------------------------
sig = 0.1;
amp = exp(-0.5*(x-0.5).^2/sig^2)';
amp0 = amp; % store initial profile

%-------------------------------------------------------------------------------
%% Set up plot for animation
%-------------------------------------------------------------------------------
% Plot the amplitude versus position, the initial amplitude
% profile, and the analytic values, with annotations
f = figure;
f.Color = 'w';
hold('on')
niceRed = [0.84,0.09,0.11];
niceBlue = [0.17,0.51,0.73];
p_Amp = plot(x,amp,'-o','color',niceRed,'LineWidth',1.5);
p_IC = plot(x,amp0,':','color',ones(1,3)*0.5,'LineWidth',1.5);
p_Ann = plot(x,amp0,'+','color',niceBlue,'LineWidth',1.5);

% Record a(x,t) matrix for visualisation
amp_xt(:,1) = amp0;
% Precompute discretization of time:
time = (0:numSteps)*tau;

%-------------------------------------------------------------------------------
% March forwards in time:
%-------------------------------------------------------------------------------
for n = 1:numSteps

    % Update the wave amplitude profile:
    amp = M*amp;

    % Calculate the profile for the exact analytic solution
    k = floor(c*time(n+1)/h); % The index the peak has reached
    xp = x([k+1:end 1:k]); % Shift x to left by k steps

    % Plot the amplitude versus position, the initial amplitude
    % profile, and the analytic values, with annotations
    p_Amp.YData = amp;
    p_Ann.XData = xp;
    p_Ann.YData = amp0;
    xlabel('Position (non-dim.)');
    ylabel('Amplitude (non-dim.)');
    title(['Time: ',num2str(time(n))]);
    drawnow();
    pause(frameUpdateLag);

    % Record a(x,t) matrix for visualisation
    amp_xt(:,n+1) = amp;
end

%-------------------------------------------------------------------------------
% Add legend
handle = legend('Numerical solution','Initial profile','Analytic solution');
set(handle,'Box','off','Location','NorthWest')

% Visualisation of amplitude versus position and time. MATLAB
% uses a surfc(x(i),y(j),z(j,i)) convention, hence the transpose
figure(2);
colormap(summer);
surfc(x,time,amp_xt');
shading interp;
xlabel('Position');
ylabel('Time');
zlabel('Amplitude');
