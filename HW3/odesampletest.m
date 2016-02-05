
% odesampletest
% test odesample for various tolerances
%
% From  http://www.amath.washington.edu/~rjl/fdmbook/chapter5  (2007)

clear all
ODE45 = ' ode45';
ODE113 = 'ode113';
global fcnevals
fprintf('Results or %s Solver', ODE113)
disp(' ')
disp('       tol          max error    f evaluations')
disp(' ')
for tol = logspace(-1,-13,13)
   %odesample(tol)
   err = Problem5_8_a(tol, 'off', ODE113);
   disp(sprintf('  %12.3e   %12.3e   %7i',tol, err,fcnevals))
   end
   
fprintf('Results or %s Solver', ODE45)
disp(' ')
disp('       tol          max error    f evaluations')
disp(' ')
for tol = logspace(-1,-13,13)
   %odesample(tol)
   err = Problem5_8_a(tol, 'off', ODE45);
   disp(sprintf('  %12.3e   %12.3e   %7i',tol, err,fcnevals))
   end
   