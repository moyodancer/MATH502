close all
clear all

%% Figure 2.4a
T = 2*pi;
alpha = .7;
beta = .7;
figure(1)
hold on
m = 20;
t = linspace(0, T, m+2);
theta(:) = (.7*cos(t) +.5*sin(t));
theta(1) = alpha;
theta(m+2) = beta;
for i = 1:5
    plot( t, theta)
    J = makeJ(theta, T);
    G = makeG_nonlinPendulum(theta, T);
    delta = J\-G';
    delta = [0; delta; 0];
    theta = theta + delta';
end
title('Figure 2.4a Theta = .7*cos(t) + .5*sin(t)')
ylabel('theta (radians)')
xlabel('time (s)')
legend('k = 1', 'k = 2', 'k = 3', 'k = 4', 'k = 5')

%% Figure 2.4b
figure(2)
hold on
theta(:) = t./t*.7;
theta(1) = alpha;
theta(m+2) = beta;
for i = 1:5
    plot( t, theta)
    J = makeJ(theta, T);
    G = makeG_nonlinPendulum(theta, T);
    delta = J\-G';
    delta = [0; delta; 0];
    theta = theta + delta';
end
title('Figure 2.4b Theta = .7')
ylabel('theta (radians)')
xlabel('time (s)')
legend('k = 1', 'k = 2', 'k = 3', 'k = 4', 'k = 5')
%% Figure 2.5
figure(6)
hold on
theta(:) = .7+sin(t/2);
theta(1) = alpha;
theta(m+2) = beta;
for i = 1:5
    plot( t, theta)
    J = makeJ(theta, T);
    G = makeG_nonlinPendulum(theta, T);
    delta = J\-G';
    delta = [0; delta; 0];
    theta = theta + delta';
end
title('Figure 2.5 Theta = .7+sin(t/2)')
ylabel('theta (radians)')
xlabel('time (s)')
legend('k = 1', 'k = 2', 'k = 3', 'k = 4', 'k = 5')
%% Other Solution
figure(3)
hold on
theta(:) = .7./t;
theta(1) = alpha;
theta(m+2) = beta;
for i = 1:5
    plot( t, theta)
    J = makeJ(theta, T);
    G = makeG_nonlinPendulum(theta, T);
    delta = J\-G';
    delta = [0; delta; 0];
    theta = theta + delta';
end
title('Alternate Solution: theta = .7/t')
ylabel('theta (radians)')
xlabel('time (s)')
legend('k = 1', 'k = 2', 'k = 3', 'k = 4', 'k = 5')

%% Large T
clear t theta m T J G delta
figure(4)
hold on
m = 100;
T = 20;
t = linspace(0, T, m+2);
theta(:) = sin(pi*t/T);
theta(1) = alpha;
theta(m+2) = beta;
for i = 1:5
    plot( t, theta)
    J = makeJ(theta, T);
    G = makeG_nonlinPendulum(theta, T);
    delta = J\-G';
    delta = [0; delta; 0];
    theta = theta + delta';
end
title('Large T: theta = sin(pi*t/T)')
ylabel('theta (radians)')
xlabel('time (s)')
legend('k = 1', 'k = 2', 'k = 3', 'k = 4', 'k = 5')
%% Max Theta
clear t theta m T J G delta
figure(5)
hold on
m = 100;
for T = 20:20:100
    t = linspace(0, T, m+2);
    theta(:) = sin(pi*t/T);
    theta(1) = alpha;
    theta(m+2) = beta;
    J = makeJ(theta, T);
    G = makeG_nonlinPendulum(theta, T);
    delta = J\-G';
    delta = [0; delta; 0];
    theta = theta + delta';
    plot( t, theta)
    clear t theta
end
title('Max Theta')
ylabel('theta (radians)')
xlabel('time (s)')
legend('T = 20','T= 40', 'T = 60', 'T = 80', 'T = 100')