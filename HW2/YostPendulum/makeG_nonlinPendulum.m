function [ G ] = makeG_nonlinPendulum( theta, T )
%Create the G matrix for the nonlinear pendulum problem
%   @param theta    Array of theta angles rad (1xm)
m = length(theta)-2; %interior points
h = T/(m+1);
for i = 2:m+1
    G(i) = (theta(i-1) -2*theta(i) + theta(i+1))/h^2 + sin(theta(i));
end
%G(1) = 1;
%G(m+2) = 1;
G = G(2:m+1);

end

