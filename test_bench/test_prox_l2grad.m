function errors = test_prox_l2grad()

errors =  0;
% Test Fourier

N = 10;
M = 20;

% 1D Test
x = zeros(N,1);
x(1) = 1;


A = @(x) difc(circshift(x,1));
At = @(x) -difc(x);
% 
% figure
% n=(0:N-1)';
% plot(n,abs(fft(A(x))),'xr',n,sqrt(2-2*cos(2*pi*n/N)),'bo')

paraml2.tight = 0;
paraml2.A = A;
paraml2.At = At;
paraml2.maxit = 100;
paraml2.nu = 4;
paraml2.verbose = 0;
paraml2.tol = 1e-10;
sol_opt = prox_l2(x,1,paraml2);

paraml2grad.verbose = 0;
sol_prox = prox_l2grad(x,1,paraml2grad);
% 
% figure
% plot(1:N,sol_opt,'xb',1:N,sol_prox,'ro')
% legend('opt','exact')

errors = errors + norm(sol_opt-sol_prox,'fro')>1e-5;

%%
% 2D Test

z = zeros(N,M);
z(1) = 1;

tau = 1;

A1 = @(x) difc(circshift(x,1));
A1t = @(x) -difc(x);
A2 = @(x) transpose(difc(circshift(transpose(x),1)));
A2t = @(x) transpose(-difc(transpose(x)));
f1.grad = @(x) (x-z) + tau*2*A1t(A1(x)) + tau*2*A2t(A2(x));
f1.eval = @(x) 0.5 * norm(x-z,'fro')^2 + tau* norm(A1(x),'fro')^2 + tau* norm(A2(x),'fro')^2;
f1.beta = tau*4+tau*4+1;

paramsolver.maxit = 1000;
paramsolver.verbose = 0;
paramsolver.tol = 1e-10;
sol_opt = solvep(z,f1,paramsolver);
paramprox.d2 = 1;
paramprox.verbose = 0;
sol_prox = prox_l2grad(z,tau,paramprox);
% 
% figure
% subplot(121)
% imagesc(sol_opt);
% colorbar
% subplot(122)
% imagesc(sol_prox);
% colorbar

errors = errors + norm(sol_opt-sol_prox,'fro')>1e-5;


% 
% n=(0:N-1)';
% m=(0:M-1);
% % b) values
% eig_n = (2-2*cos(2*pi*n/N));
% eig_m = (2-2*cos(2*pi*m/M));
% % c) Compute the radius for the kernel
% rho = repmat(eig_n,1,M) + repmat(eig_m,N,1);
% h=@(t) 1./(1+2*tau*t);
% 
% figure
% subplot(121)
% imagesc(abs(fft2(A1(x))))
% colorbar
% 
% subplot(122)
% imagesc(sqrt(repmat(eig_n,1,M)))
% colorbar
% 
% figure
% subplot(121)
% imagesc(abs(fft2(sol_opt)))
% colorbar
% subplot(122)
% imagesc(h(rho))
% colorbar

if errors
    warning('Errors in test_prox_l2grad')
end

end