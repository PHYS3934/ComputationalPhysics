% exp_euler.m
% Integrate dx/dt = x using Euler's method

% Clear the workspace and show only a few digits
clear('all');
format('short');

% Number of steps and time step
numSteps = 20;
tau = 1/numSteps;

% Initial values
t = 0;
x = 1;

% Show values of independent and dependent variables
fprintf('\n      t       x\n');
fprintf('%7.5g %7.5g\n',t,x);

%-------------------------------------------------------------------------------
% Euler's method integration
for n = 1:numSteps

    % One step of Euler:
    % rhs_exp(x) = x.
    f = rhs_exp(x);
    x = x + tau*f;
    t = t + tau;

    % Show values of independent and dependent variables:
    fprintf('%7.5g %7.5g\n',t,x);
end

%-------------------------------------------------------------------------------
% Display percentage error:
fprintf(1,'Error: %g%%\n',100*abs(x-exp(1))/exp(1));
