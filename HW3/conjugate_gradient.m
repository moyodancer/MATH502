%function u = conjugate_gradient(A,f,tol)
%
%   Example:  
%   x = conjugate_gradient(A,b,tol) 
% 
f = F;
tol = 1e-5;
MAXITS = length(f);

u = 0*f;
r = f-A*u;
p = r;
for k = 1:MAXITS
    w = A*p;
    alpha = (r'*r)/(p'*w);
    unew = u+alpha*p;
    rnew = r - alpha*w;
    if( norm(rnew) < tol ),
        fprintf('Converged! its= %7.0f, tol=%10.3e\n', [k tol]);
        return;
    end
    beta = (rnew'*rnew)/(r'*r);
    p = rnew + beta*p;
    r = rnew;
    u = unew;
end
fprintf('Caution: CG went to max iterations without converging!\n');
fprintf('MAXITS = %7.0f, tol =%10.3e\n', [MAXITS tol]);

%end
