close all
clear all

%% Problem 2
% Morgan Yost
% Math 502 Take Home Midterm
%%
% Set Up
T = 5;
alpha = 0;
beta = 2*pi;
force = .19;
m = 100;
t = linspace(0, T, m+2);
theta(:) = (sin(t));
theta(1) = alpha;
theta(m+2) = beta;
err = 10;
deltaOld = 1;
tol = 1e-6;
count = 0;
while(err >tol)
    J = makeJ(theta, T);
    G = makeG_nonlinPendulum(theta, T, force);
    delta = J\-G';
    delta = [0; delta; 0];
    theta = theta + delta';
    err = norm(delta-deltaOld);
    deltaOld = delta;
    count = count+1;

end
%fprintf('Completed in %d iterations\n', count');
figure(1)
plot(t, theta)
title('Figure 2.4a Guess \theta = sin(t)')
ylabel('theta (radians)')
xlabel('time (s)')

%Get the initial angulat velocity
stencil = [ 0 1];
c = fdcoeffF(1,0,stencil);
h = T/(m+1);
x0 = 0;
theta_primeInit = 1/h*sum(c.*theta(x0+stencil+1));
fprintf('Initial Angular Velocity \n %f rad/s\n', theta_primeInit)

%Sanity Check: Just want to look at the mean
stencil = [-1 0 1]; %one sided second order
c = fdcoeffF(1,0,stencil);
i = 1;
for x0 = 2:length(theta)-2
    theta_prime(i) = 1/h*sum(c.*theta(x0+stencil));
    i = i+1;
end

%Get the final angular velocity
stencil = [-1 0]; %one sided second order
c = fdcoeffF(1,0,stencil);
x0 = length(theta);
theta_primeFinal = 1/h*sum(c.*theta(x0+stencil));
fprintf('Final Angular Velocity \n %f rad/s\n', theta_primeFinal)


%% Pretty Animation of the ride
% tic
% for i = 1:length(theta)
%     pointX = [0 cos(theta(i))];
%     pointY = [0 sin(theta(i))];
%     clf 
%     plot(pointX, pointY )
%     axis([-1 1 -1 1])
%   drawnow
%   pause(1/(m+1))
% end
% toc
