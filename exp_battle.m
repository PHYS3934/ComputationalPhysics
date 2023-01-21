% exp_battle.m
% Integrate dx/dt = x
% Compare Euler and RK4

% Clear the workspace
clear('all');

% Number of steps and time step
numSteps = 5;
tau = 1/numSteps;

% Initial values
t = zeros(numSteps+1,1);
x = zeros(numSteps+1,2);
x(1,:) = 1;

% Show values of independent and dependent variables
% fprintf('\n      t       Euler x        RK4 x\n');
% fprintf('%7.5g %7.5g  %7.5g\n',t(1),x(1,1),x(1,2));

f = figure('color','w');
hold('on')
xlabel('t')
ylabel('x')
h_Euler = plot(t(1),x(1,1),'o-r')
h_RK4 = plot(t(1),x(1,2),'o-b')

%-------------------------------------------------------------------------------
% Euler's method integration
for n = 1:numSteps

    % One step of Euler:
    f = rhs_exp(x(n,1));
    x(n+1,1) = x(n,1) + tau*f;

    % One step of RK4:
    f1 = rhs_exp(x(n,2));
    f2 = rhs_exp(x(n,2) + 0.5*tau*f1);
    f3 = rhs_exp(x(n,2) + 0.5*tau*f2);
    f4 = rhs_exp(x(n,2) + tau*f3);

    x(n+1,2) = x(n,2) + tau*(f1 + 2*f2 + 2*f3 + f4)/6;

    % Increment time:
    t(n+1) = t(n) + tau;

    % Print out values of independent and dependent variables:
    tRange = linspace(0,t(n+1),100);
    plot(tRange,exp(tRange),'-k')
    h_Euler.XData = t(1:n+1);
    h_Euler.YData = x(1:n+1,1);
    h_RK4.XData = t(1:n+1);
    h_RK4.YData = x(1:n+1,2);
    xlim([0,1])
    ylim([1,exp(1)])
    drawnow()
    input('')

end
