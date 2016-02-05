clear all
addpath ../fdmbook
%% 2.2: Green's Function with Neumann and Drichlet Boundary Conditions
% Equation 2.54
h = .25;
A = 1/h^2*[ -h h 0 0 0; 1 -2 1 0 0; 0 1 -2 1 0; 0 0 1 -2 1; 0 0 0 0 h^2];

% Green's function to get A inverse
m = 3;
x = 0:1/(m+1):1;

for col = 0:m+1
    Ainv(:, col+1) = getG2_2(col, m, h, x);
end
Ainv

%% 2.3: Solvability Condition for Nuemann Problem
clear all
% Equation 2.58
h = .25;
A =  1/h^2*[ -h h 0 0 0; 1 -2 1 0 0; 0 1 -2 1 0; 0 0 1 -2 1; 0 0 0 h -h];
x = null(A')
% Equation 2.62
f = sum((A*x).*[h/2 h h h h/2]') %??

%% 2.4: Modifying BVP code
clear all
close all
% Part A
bvp_2_2_4
% Part B
bvp_4_2_4

%% 2.6 Ill Posed BVP
clear all 
close all
ax = 0;
bx = 1;
alpha = 2;   % Dirichlet boundary condition at ax
beta = 3;     % Dirichlet boundary condtion at bx
bvp_2_2_6
%%
clear all
ax = 0;
bx = pi;
count = 1;
check = [0 pi/4 pi/2 3*pi/4 pi];
x = linspace(0, pi);
syms A B
for alpha = check
    for beta = check
        X = bvp_check(ax, bx, alpha, beta, A, B);
        if X(1)~= 0 && X(2) ~= 0
            keep(:, count) = [alpha; beta];
            count = count +1;
            sln(:, count) = X(1)*cos(x) + X(2)*sin(x);
        end
    end
end
figure (1)
hold on
for i = 1:count
    plot(x, sln(:, i))
end
title('Solution Space for a = 0 and b = pi')
xlabel('x')
ylabel('u(x)')
%%
clear all
ax = 0;
bx = pi;
alpha = 1;   % Dirichlet boundary condition at ax
beta = -1;     % Dirichlet boundary condtion at bx
bvp_2_2_6
%%
clear all
ax = 0;
bx = pi;
alpha = 1;   % Dirichlet boundary condition at ax
beta = 1;     % Dirichlet boundary condtion at bx
bvp_2_2_6

