function [ errors ] = test_ppxa_and_other( )
%TEST_PROX_L2 This function test the function prox_l2
errors=0;
tol=10e-10;
errors=errors+test_proj(tol);
errors=errors+test_prox_l21(tol);
errors=errors+test_svst_l21(tol);
errors=errors+test_ppxa_only(10e-5);

end

function [errors]= test_svst_l21(tol)
    A=rand(100);
    A=A*A;
    lambda=0.1;
    param.verbose=0;
    [U, Sigma, V, ~, ~] = svst(A, lambda, param);
    
    p2 = U * Sigma * V';

    p3=prox_nuclearnorm(A,lambda,param);

    
    if norm(p3-p2)/norm(p3)<tol
        fprintf('  Test svst OK\n')
        errors=0;
    else
        fprintf('  Test svst Pas OK!!!!!!!!!!!!!!!!\n')
        
        errors=1;
    end
end

function [errors]= test_prox_l21(tol)
    
    %****Prox L2/L1********************************************
    lambda=0.1;
    N=100;
    x=rand(N,1);
    n2=10;
    param_l21.verbose=0;
    param_l21.g_d=(reshape(reshape(1:N,N/n2,n2)',N,1))';
    param_l21.g_t=(N/n2*ones(n2,1))';
    
    p2 = reshape(x,[], n2);
    l = sqrt(sum(p2.^2 , 2))+ eps;
    l = max(l - lambda,0) ./ l;
    p2 = p2 .* repmat(l,1,n2);
    
    p2 = p2(:);
    p3=prox_l21(x,lambda,param_l21);
   
    
    if norm(p3-p2)<tol
        fprintf('  Test prox_l21 OK\n')
        errors=0;
    else
        fprintf('  Test prox_l21 Pas OK!!!!!!!!!!!!!!!!\n')
        errors=1;
        norm(p3-p2)/norm(p3)
    end
    
end

function [errors]= test_proj(tol)
    or=rand(100,1);
    
    Am=randn(500,100);
    A = @(x) Am * x;
    At = @(x) Am' * x;
    
    y=A(or);
    
    paramL2.epsilon = 0.001;
    paramL2.max_iter = 200;
    paramL2.verbose = 0;
    paramL2.tight = 0;
    paramL2.nu = norm(Am)^2;  
    
    xin=zeros(size(or));
    
    p1 = proj_b2_test(xin, y, A, At, paramL2);
    
    paramL2.A=A;
    paramL2.At=At;
    paramL2.y=y;
    paramL2.method='ISTA';
    
    p2= proj_b2(xin,0,paramL2);
    
    if norm(p1-p2)<tol
        fprintf('  Test proj OK\n')
        errors=0;
    else
        fprintf('  Test proj Pas OK!!!!!!!!!!!!!!!!\n')
        errors=1;
    end
    
    
end

function [errors]= test_ppxa_only(tol)
    n1 = 10;
    n2 = 10;
    r = 2;   % rank of data matrix
    k1 = 5; % Joint-sparsity level
    m =  4*r*(n1+n2-r) ; % nb_measurements = 4 * (data degrees of freedom) 


    %--- Random measurement matrix--------
    A_mtx = 1/sqrt(m) * randn(m, n1*n2);
    A = @(x) A_mtx * x;
    At = @(x) A_mtx' * x;

    % joint-sparse lowrank data generation
    x = zeros(n1,n2);
    nz_ind = randsample (n1, k1);
    x(nz_ind,:) = randn(k1,r)*randn(r,n2);

    %-- Compressed sampling-----
    y = A(x(:));
    z = 0;         %1e-8*randn(nb_meas*J,1); % Additive white gaussian noise.
    y = y + z;

    % ****** Recovery of from CS measurements.
    %fprintf('Recovery commenced...\n')
    param.alpha =  sqrt(2*r/k1);             % Regularization parameter
    param.epsilon = 0;                       % For noisy data (Fidelity bound)
    %param.lambda = .1 * svds(x(nz_ind,:),1); % Parameter of singular value soft thresholding.
    param.tight = 0;                         % (0 if A is not tight-frame) (1 if A is a tight-frame)
    param.nu = norm(A_mtx)^2;                % Upper bound on the spectral norm of the forward operator 
    param.verbose = 0;                       % Graphical feedback at each main iteration.
    param.tol = 0;
    param.MaxIter = 50;
    param.maxit=50;
    param.r_init=n1*n2;
    xhat = LR_plus_MixNorm_PPXA(y, n1, A, At, param); % LR-JS coefficiennts
    xhat = reshape(xhat, n1, n2);
    
    % projection
    paramL2.epsilon = param.epsilon;
    paramL2.max_iter = 200;
    paramL2.verbose =0;
    paramL2.tight = param.tight;
    paramL2.nu = param.nu;
    paramL2.verbose=0;
    
    F{1}.prox=@(x,T) proj_b2_test(x, y, A, At, paramL2);
    F{1}.eval=@(x) eps;
    
    %prox l21
    param_l21.g_d=(reshape(reshape(1:n1*n2,n1,n2)',n1*n2,1))';
    param_l21.g_t=(n1*ones(n2,1))';
    param_l21.verbose=0;
    
    F{2}.prox=@(x,T) (prox_l21(x,param.alpha*T,param_l21));
    F{2}.eval=@(x) param.alpha*norm_l21(x,param_l21.g_d,param_l21.g_t);
    
    %prox nuclearnorm
    param_svt.verbose=0;
    F{3}.prox=@(x,T) reshape(prox_nuclearnorm(reshape(x,n1,n2),T,param_svt),n1*n2,1);
    F{3}.eval=@(x) sum((svd(reshape(x,n1,n2))));
    
    x0 = zeros(n1*n2,1);
  
    xhat2=ppxa(x0,F,param);
    xhat2 = reshape(xhat2, n1, n2);

    if norm(xhat2-xhat)/norm(xhat)<tol
        errors=0;
        fprintf('  Test ppxa OK\n')
    else
        errors=1;
        norm(xhat2-xhat)/norm(xhat)
        fprintf('  Test ppxa pas OK!!!!!!!!\n')
    end

end


function [x] = LR_plus_MixNorm_PPXA(y, n1, A, At, param)

% This program solves the follwoing optimization problem,
%
% argmin_x  |x|_* + alpha |x|_{2,1}
%  s.t.
% |y - A(x)|_F <= epsilon

%
% 1- We use PPXA Algorithm. 
% Ref: P. L. Combettes and J. C. Pesquet, ?Proximal splitting methods in
% signal processing, 2011.
%
% 2- We use SVT package for singular value soft thresholding:
% http://www-stat.stanford.edu/~candes/svt/
% 4- A forward operator, At its tarnspose.
% 5- param.tight = 1 (if A is a tight frame), 0 (else).
% 6- param.nu = norm(A)^2;  The spectral norm of A
% 7- param.epsilon = epsilon
%
% Please cite our related papers to this algorithm:
% 
% M. Golbabaee and P. Vandergheynst. "Guaranteed recovery of a low-rank and
% joint-sparse matrix from incomplete and noisy measurements,"
% Workshop on Signal Processing with Adaptive Sparse Structured Representations (SPARS11), 
% Edinburgh, UK, June 2011.
%
%
% M. Golbabaee and P. Vandergheynst. "HYPERSPECTRAL IMAGE COMPRESSED SENSING
% VIA LOW-RANK AND JOINT-SPARSE MATRIX RECOVERY," 
% Submitted to ICASSP 2012.
%
% Author: Mohammad Golbabaee
% e-mail: mohammad.golbabaei@epfl.ch
% 08.06.2011, Lausanne, Switzerland

if ~isfield(param, 'MaxIter'), param.MaxIter = 200; end
if ~isfield(param, 'tol'), param.tol = 1e-6; end
if ~isfield(param, 'lambda'), param.lambda = 0.99; end
if ~isfield(param, 'verbose'), param.verbose = 1; end

paramL2.epsilon = param.epsilon;
paramL2.max_iter = 200;
paramL2.verbose = 0;
paramL2.tight = param.tight;
paramL2.nu = param.nu;
alpha = param.alpha;

[n2] = size(At(y),1) / n1;
if ~isfield(param, 'r_init'), param.r_init = min(n1,n2)-1; end
%***Initialization************
x = zeros(n1*n2,1);
z1 = x;
z2 = x;
z3 = x;
Sigma = [];
svst_iter = [];
r = param.r_init;
Iter = 0;
x_rel_err = 1.1*param.tol+eps;
%**** Main Douglas-Rachford Iteration Loop *************************************************
while ( (Iter < param.MaxIter) && (x_rel_err > param.tol))
    
    Iter=Iter+1;
    
    %****Prox L2 ********************************************
    p1 = proj_b2_test(z1, y, A, At, paramL2);
    %****Prox L2/L1********************************************
    
    p2 = reshape(z2,[], n2);
    l = sqrt(sum(p2.^2 , 2))+ eps;
    l = max(l - param.lambda*(alpha),0) ./ l;
    p2 = p2 .* repmat(l,1,n2);
    p2 = p2(:);
    
    %*** Prox_lambda/*f1 (Singular Value Soft Thesholding)***
    if Iter > 1
        z3 = reshape(z3, [], n2);
        param_svst.s_initial = min ( [r + 1 , n1 , n2] );
        [U, Sigma, V, r, svst_iter] = svst(z3, param.lambda, param_svst);
        p3 = U * Sigma * V';
        z3 = z3(:);
        p3 = p3(:);
    elseif Iter == 1
        p3 = x;
    end
    
    
    %********AVRAGING***************************************
    x_tmp = (p1 + p2 + p3)/3;
    
    
    %*** Global Error update *****************************************
    %x_rel_err = norm(x-x_tmp)/norm(x);
    
    
    x_tmp = reshape(x_tmp,n1,n2);
    obj_tmp = sum(svds(x_tmp,r)) + alpha*sum(sqrt(sum(x_tmp,2).^2));
    x = reshape(x,n1,n2);
    obj = sum(svds(x,r)) + alpha*sum(sqrt(sum(x.^2,2)));
    x_tmp = x_tmp(:);
    x = x(:);    
    x_rel_err = abs(obj-obj_tmp)/obj_tmp;
    
    %x = x_tmp;
    
    if param.verbose >= 1
        fprintf('Main_Iter: %i, Rank= %i, |x|_* = %e, Maxiter_svst=%i , Sol_update=%e , curr_norm=%e \n',...
            Iter, r, sum(diag(Sigma)), svst_iter , x_rel_err, obj);
    end
    
    %***** Verifying theGlobal stopping criterion***********************
    if (x_rel_err < param.tol)
        %crit_BPDN = 'REL_NORM';
        break;
    elseif (Iter > param.MaxIter)
        
        break
    end
    
    %*** Solution Updates***********************************************
    z1 = z1 + 2*x_tmp - x -p1;
    z2 = z2 + 2*x_tmp - x -p2;
    z3 = z3 + 2*x_tmp - x -p3;
    
    x = x_tmp;
    
end
end
