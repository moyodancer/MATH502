function [ J ] = makeJ( theta, T )
%Create the Jacobian Matrix to solve the nonlinear pendulum problem using
%Newton's method.
%   @param theta    Array of theta angles rad (1xm)

m = length(theta)-2; %interior points
h = T/(m+1);
for i = 2:m+1
    J(i, i-1) = 1/h^2;
    J(i, i) = (-2/h^2 +cos(theta(i)));
    J(i, i+1) = 1/h^2;
end
% J(1, 1) = 1/h^2*(-2+(h^2)*cos(theta(1)));
% J(1, 2) = 1/h^2;
% J(m+2, m+1) = 1/h^2;
% J(m+2, m+2) = 1/h^2*(-2+(h^2)*cos(theta(m)));
J = J(2:m+1, 2:m+1);
end

