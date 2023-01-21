% diffusion_ftcs_n.m
% Solve the 1-D diffusion equation for an initial spike profile 
% with Neumann conditions using FTCS, in a matrix formulation

% Clear memory and show only a few digits
clear all; format short;

% Thermal conductivity
kappa=1;

% Time step and spatial step
tau=1e-4;
h=0.02;

% Number of time steps
nsteps=2000;

% Calculate the ratio tau/(th/2), where th is the approximate diffusion
% time for one spatial step h, and display
th=h^2/kappa;
disp(['Ratio tau/(0.5*th): ',num2str(tau/(0.5*th))]);

% Vector of x values
x=0:h:1;
L=length(x);

% Construct the matrix D associated with the second spatial 
% derivative and the boundary conditions
D=-2*eye(L);
D=D+diag(ones(L-1,1),+1)+diag(ones(L-1,1),-1);
D=kappa*tau*D/h^2;

% Construct the update matrix
A=eye(L)+D;

% Neumann boundary conditions
A(1,:)=A(2,:);
A(L,:)=A(L-1,:);

b=zeros(L,1);
b(1,1)=1;
b(L,1)=-1;

% Initial conditions
temp1=x+2;
temp1(10:L-10)=0;
temp=temp1';

% Record T(x,t) matrix for visualisation
temp_xt(:,1)=temp;

% March forwards in time
time(1)=0;
figure(1);
skip=50;
for n=1:nsteps

    % Update the time and the temperature profile
    time(n+1)=time(n)+tau;
    temp=A*temp-h*b;
    
    % Plot the temperature versus position, the initial temperature
    % profile, and the analytic values, with annotations
    if rem(n,skip) == 0
        plot(x,temp,'ro-',x,temp1,'g-*')
        xlabel('Position (non-dim.)');
        ylabel('Temperature (non-dim.)');
        handle=legend('Numerical solution','Initial profile');
        set (handle,'Box','off','Location','NorthWest')
        title(['Time: ',num2str( time(n))]);
        drawnow; pause(0.01);
    end
    
    % Record T(x,t) matrix for visualisation
    temp_xt(:,n+1)=temp;
    
end

% Visualisation of temperature versus position and time. MATLAB
% uses a surfc(x(i),y(j),z(j,i)) convention, hence the transpose
figure(2);
contourf(x,time,temp_xt');
%shading interp;
xlabel('Position'); ylabel('Time'); zlabel('Temperature');
