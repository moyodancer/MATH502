function [ G ] = makeG_nonlinPendulum( theta, T, force )
%Create the G matrix for the nonlinear pendulum problem
%   @param theta    Array of theta angles rad (1xm)
%   @param T        Period of the motion
%   @param force    Forcing constant
%   @return G       G matrix for J*delta = -G
% AUTHOR: Morgan Yost
m = length(theta)-2; %interior points
h = T/(m+1);
for i = 2:m+1
    G(i) = (theta(i-1) -2*theta(i) + theta(i+1))/h^2 + sin(theta(i))-force;
end

G = G(2:m+1);

end

