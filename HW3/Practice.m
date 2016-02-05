clear all
close all
count = 6;
count = count+1;
mx = 8;
my = 9;
ax = 0;  
bx = 1;
ay = 0;
by = 2;
hx = (bx-ax)/(mx+1);
hy = (by-ay)/(my+1);

x = linspace(ax, bx, mx+2);  % grid points x including boundaries
y = linspace(ay, by, my+2);   % grid points y including boundaries


[X,Y] = meshgrid(x,y);      % 2d arrays of x,y values
X = X';                     % transpose so that X(i,j),Y(i,j) are
Y = Y';                     % coordinates of (i,j) point

Iint = 2:mx+1;              % indices of interior points in x
Jint = 2:my+1;              % indices of interior points in y
Xint = X(Iint,Jint);       % interior points
Yint = Y(Iint,Jint);

f = @(x,y) 1.25*exp(x+y/2);         % f(x,y) function

rhs = f(Xint,Yint);  % evaluate f at interior points for right hand side
                           % rhs is modified below for boundary conditions.

utrue = exp(X+Y/2);        % true solution for test problem

% set boundary conditions around edges of usoln array:

usoln = utrue;              % use true solution for this test problem
                            % This sets full array, but only boundary values
                            % are used below.  For a problem where utrue
                            % is not known, would have to set each edge of
                            % usoln to the desired Dirichlet boundary values.


% adjust the rhs to include boundary terms:
rhs(:,1) = rhs(:,1) - usoln(Iint,1)/hy^2;
rhs(:,my) = rhs(:,my) - usoln(Iint,my+2)/hy^2;
rhs(1,:) = rhs(1,:) - usoln(1,Jint)/hx^2;
rhs(mx,:) = rhs(mx,:) - usoln(mx+2,Jint)/hx^2;


% convert the 2d grid function rhs into a column vector for rhs of system:
F = reshape(rhs,mx*my,1);

% form matrix A:
Ix = speye(mx);
Iy = speye(my);
e = ones(my,1);
Tx = spdiags([e -2*e e],[-1 0 1],mx,mx);
Ty = spdiags([0*e -2*e 0*e],[-1 0 1],mx,mx);
S = spdiags([e e],[-1 1],my,my);
A = (kron(Iy,Tx)/hx^2 + kron(Iy,Ty)/hy^2 + kron(S,Ix)/hy^2) ;


% Solve the linear system:
uvec = A\F;  

% reshape vector solution uvec as a grid function and 
% insert this interior solution into usoln for plotting purposes:
% (recall boundary conditions in usoln are already set) 

usoln(Iint,Jint) = reshape(uvec,mx,my);

% assuming true solution is known and stored in utrue:
err = max(max(abs(usoln-utrue)));   
fprintf('grid size: %dx%d\n', mx, my);
fprintf('Error relative to true solution of PDE = %10.5e \n',err)

% plot results:

figure(count)
hold on

% plot grid:
 plot(X,Y,'g');  plot(X',Y','g')

% plot solution:
contour(X,Y,usoln,30,'k')

axis([ax bx ay by])
daspect([1 1 1])
name = sprintf('Contour plot of computed solution for %dx%d rectangular grid', mx, my);
title(name)
hold off