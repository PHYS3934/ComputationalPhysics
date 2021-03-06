% laplace_analytic.m
%-------------------------------------------------------------------------------
% Evaluate the analytic solution to the Dirichlet problem for the
% Laplace equation given in Garcia
%-------------------------------------------------------------------------------

% Clear memory and show only a few digits
clear('all');
format('short');

% Spatial step and number of terms in the sum
h = 0.025;
numTerms = 25;

% Vectors of x and y values
x = 0:h:1;
y = x;

% Array of dependent variable values
L = length(x);
phi = zeros(L);
Ex = zeros(L);
Ey = zeros(L);

% Calculate the sum in the analytic solution at each grid point
nvals = 1:2:2*numTerms + 1;
for j = 1:L
    for l = 1:L
        % Terms in the numerators and denominators of the sums
        sterm = sin(nvals*pi*x(j));
        shterm = sinh(nvals*pi*y(l));
        chterm = cosh(nvals*pi*y(l));
        cterm = cos(nvals*pi*x(j));
        dterm = sinh(nvals*pi);
        phi(j,l) = 4*sum(sterm.*shterm./(nvals*pi.*dterm));
        Ex(j,l) = -4*sum(cterm.*shterm./dterm);
        Ey(j,l) = -4*sum(sterm.*chterm./dterm);
    end
end

%-------------------------------------------------------------------------------
% Plot a surface plot of phi
f1 = figure(1);
f1.Color = 'w';
set_color_map()
surfc(x,y,phi');
axis([0 1 0 1 0 1.1]);
shading('faceted');
xlabel('x');
ylabel('y');
zlabel('Potential, \phi');

%-------------------------------------------------------------------------------
% Plot the electric field
f2 = figure(2);
f2.Color = 'w';
hold('on');
set_color_map()
contourf(x,y,phi');
scale = 10;
quiver(x,y(1:L-1),Ex(:,1:L-1)',Ey(:,1:L-1)',scale,'w');
xlabel('x');
ylabel('y');
title('Electric field {\bf E} and contours of \phi');
axis([0 1 0 1]);
hold('off');
