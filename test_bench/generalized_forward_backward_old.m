function [sol, info,objective] = generalized_forward_backward_old(x_0, F, f , param)
% Start the time counter
t1 = tic;

% Optional input arguments
if nargin<4, param=struct; end

if ~isfield(param, 'weights'), param.weights=ones(size(F,2),1); end
if ~isfield(param, 'verbose'), param.verbose=1 ; end
if ~isfield(param, 'maxit'), param.maxit=100 ; end
if ~isfield(param, 'tol'), param.tol=1e-3 ; end
if ~isfield(param, 'gamma'), param.gamma=1 ; end
if ~isfield(param, 'lambda'), param.lambda=0.99 ; end

if nargin<3 
    f.grad=@(x) 2*x;
    f.eval=@(x) norm(x(:)-x_0(:),2)^2;  
end


% Normalizing the weights
param.weights= param.weights/sum(param.weights);

% Number of function
l=size(F,2);



% test the evaluate function
f = test_eval(f);
F = test_eval(F);


% Definition of the gradient function
grad = @(y) f.grad(y);

% Algorithm - Initialisation
z=zeros([l,size(x_0)]);


for ii=1:l
    z(ii,:,:,:,:)=x_0;
end


sol=x_0;
curr_norm = f.eval(sol)+norm_sumg(sol,F);
[~,~,prev_norm,iter,objective,~] = convergence_test_old(curr_norm);

% Algorithm - Loop

while 1
    
    %
    if param.verbose > 1
        fprintf('Iteration %i:\n', iter);
    end
    temp_grad=grad(sol);
    for ii=1:length(F)
       temp=2*sol-reshape(z(ii,:,:,:,:),size(x_0))-param.lambda*temp_grad;
       z(ii,:,:,:,:) = reshape(z(ii,:,:,:,:),size(x_0))+ param.gamma*(F{ii}.prox(temp,1/param.weights(ii)*param.lambda)-sol);
    end
    
    sol=zeros(size(x_0));
    for ii=1:l
        sol=sol+param.weights(ii)* reshape(z(ii,:,:,:,:),size(x_0));
    end
    
    % Global stopping criterion

    curr_norm = f.eval(sol)+norm_sumg(sol,F);
    [stop,rel_norm,prev_norm,iter,objective,crit] = convergence_test_old(curr_norm,prev_norm,iter,objective,param);
    [sol, param] = post_process_old(sol, iter, curr_norm, prev_norm, objective, param);
    if stop
        break;
    end
    if param.verbose > 1
        fprintf(' Norm of the general objectiv function: ||f|| = %e, rel_norm = %e\n', ...
            curr_norm, rel_norm);
    end


  
    
end





% Calculation of the norm
norm_G=f.eval(sol)+norm_sumg(sol,F);

% Log after the calculous of the prox
if param.verbose > 1
    fprintf('  Generalized forward backward: Sum_k ||G_k(x)|| = %e\n', norm_G);


    % Stopping criterion
    fprintf(' %i iterations\n', iter);
    fprintf(' Stopping criterion: %s \n\n', crit);

elseif param.verbose==1
    fprintf('  Generalized forward backward: Sum_k ||G_k(x)|| = %e', norm_G);


    % Stopping criterion
    fprintf(', %i iter', iter);
    fprintf(', crit: %s \n', crit);
end

info.algo=mfilename;
info.iter=iter;
info.final_eval=curr_norm;
info.crit=crit;
info.time=toc(t1);
info.rel_norm=rel_norm;

end