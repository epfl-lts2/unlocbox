function [sol,info] = prox_linf1(x, gamma , param)
%PROX_LINF1 Proximal operator with L1inf norm
%   Usage:  sol = prox_linf1(x, gamma, param)
%           [sol,info] = prox_linf1(x, gamma, param)
%
%   Input parameters:
%         x     : Input signal.
%         gamma : Regularization parameter.
%         param : Structure of parameters (optional)
%   Output parameters:
%         sol   : Solution.
%         info  : Structure summarizing informations at convergence
%
%   `prox_Linf1(x, gamma, param)` solves:
%
%   .. sol = argmin_{z} 0.5*||x - z||_2^2 + gamma * ||A x||_inf1
%
%   .. math::  sol = \min_{z} \frac{1}{2} \|x - z\|_2^2 + \gamma  \|A x\|_{\infty 1}
%
%   The easiest way to use this proximla operator is to give a matrix x as
%   imput. In this case, the sup norm will be computed over the lines (2nd
%   dimention) and the one norm will be computed over the rows (1st
%   dimention).
%
%   *param* is a Matlab structure containing the following fields: 
%
%   * *param.weights1* : weights for a weighted L1inf-norm works on the
%     norm L1 (default = 1) (Experimental)
%
%   * *param.weights2* : weights for a weighted L1inf-norm works on the sup
%     nom (default = 1)  (Experimental)
%
%   * *param.g_d*, *param.g_t* are the group vectors. If you give a matrix,
%     do not set those parameters.
%
%     *param.g_d* contains the indices of the elements to be grouped and
%     *param.g_t* the size of the different groups. 
%
%     Warning: *param.g_d* and *param.g_t* have to be row vector!     
%     
%     Example: suppose x=[x1 x2 x3 x4 x5 x6] 
%                  and Group 1: [x1 x2 x4 x5] 
%                      group 2: [x3 x6]
%              
%     In matlab:: 
%
%           param.g_d = [1 2 4 5 3 6]; param.g_t=[4 2];
%
%     Also this is also possible::
%
%           param.g_d = [4 5 3 6 1 2]; param.g_t=[2 4]; 
%
%   * *param.multi_group*: in order to group component in a not disjoint
%     manner, it is possible to use the multi_group option.
%     *param.multi_group* is now set automatically by the function. 
%
%     Overlaping group:
%     In order to make overlapping group just give a vector of g_d, g_b
%     and g_t. Example::
%       
%           param.g_d=[g_d1; g_d2; ...; g_dn];
%           param.g_t=[g_t1; g_t2; ...; g_tn];
%
%     Warning! There must be no overlap in *g_d1*, *g_d2*,... *g_dn*
%
%   * *param.verbose* : 0 no log, 1 a summary at convergence, 2 print main
%     steps (default: 1)
%
%
%   info is a Matlab structure containing the following fields:
%
%   * *infos.algo* : Algorithm used
%
%   * *info.iter* : Number of iteration
%
%   * *info.time* : Time of exectution of the function in sec.
%
%   * *info.final_eval* : Final evaluation of the function
%
%   * *info.crit* : Stopping critterion used 
%
%
%   See also:  prox_l1 prox_l21 proj_b1 prox_sumg
%
%   Demos: demo_compress_sensing4
%
%   References: bach2011optimization


% Author: Nathanael Perraudin
% Date: November 2012
% Testing: test_mixed_sparsity

% Start the time counter
t1 = tic;

% Reshape x if not a row vector
t=size(x);

% Optional input arguments
if nargin<3, param=struct; end

if ~isfield(param, 'g_d'),    param.g_d = 1:numel(x); end
if ~isfield(param, 'g_t')
    if numel(x) == size(x,1)*size(x,2); % matrix case
        param.g_t = size(x,2)*ones(1,size(x,1)); 
    else
        param.g_t = ones(1,numel(x)); 
    end
end
if ~isfield(param, 'weights1'), param.weights1=ones(length(param.g_t),1) ; end
if ~isfield(param, 'weights2'), param.weights2=ones(numel(x),1) ; end
if ~isfield(param, 'verbose'), param.verbose=1 ; end
if ~isfield(param, 'multi_group'), param.multi_group=0 ; end

% Test of gamma
if test_gamma(gamma)
    sol = x;
    info.algo=mfilename;
    info.iter=0;
    info.final_eval=0;
    info.crit='--';
    info.time=toc(t1);
    return; 
end

% Test of the weights
param.weights1=test_weights(param.weights1);
param.weights2=test_weights(param.weights2);

test_multigroup(x,param.g_d,param.g_t);

% test if there are more than one group
if size(param.g_d,1)>1
    param.multi_group=1;
end



if param.multi_group==0;
    % Number of group
    l=length(param.g_t);
   
    % Test if all the group have the same size
    if max(param.g_t)==min(param.g_t),
        
        
        % reshape x in a useful manner
        X=transpose(x);
        X=X(param.g_d);
        X=transpose(reshape(X,numel(x)/l,l));
        W2=transpose(reshape(param.weights2(param.g_d),numel(x)/l,l));
        
        % soft thresholding
        temp=W2.*X;
        S=gamma./max(abs(temp),[],2).*param.weights1;
        
        sol=soft_thresholdb( temp,repmat(S,1,max(param.g_t)));
        
        % handle size for row vector
        if size(x,1)*size(x,2) == length(x)
            sol = transpose(sol);
        end
              
        %reconstruct the solution
        sol(param.g_d)=sol;
        
        
        

        
    else % group of different size
        sol=zeros(size(x));
        indice=0;
        
        xp=x;
        xp=xp(param.g_d);
        W2=param.weights2;
        W2=W2(param.g_d);
        for i=1:l
           temp=xp(indice+1:indice+param.g_t(i));
           w2=W2(indice+1:indice+param.g_t(i));
           s=soft_thresholdb( temp, ...
               gamma/norm(w2.*temp,Inf)*param.weights1(i));
           sol(indice+1:indice+param.g_t(i))=s;
           indice=indice+param.g_t(i);
           
        end
        %reconstruct the solution
        sol(param.g_d)=sol;
    end
    
    
    
    norm_Linf1=norm_linf1(x,param.g_d,param.g_t,param.weights2,param.weights1);
    % Log after the calculous of the prox
    if param.verbose >= 1
        fprintf('  prox_L12: ||x||_12 = %e\n', norm_Linf1);
    end

else
    r=size(param.g_t,1);
    
    % Parameter for the prox

    G=[];
    
    for k=1:r
        param3.g_t=param.g_t(k,:);
        param3.g_d=param.g_d(k,:);
        param3.multi_group=0;
        param3.verbose=0;
        g.prox=@(x,T) prox_linf1(x,T,param3);
        g.eval=@(x) norm_linf1(x,param.g_d(k,:),param.g_t(k,:));
        G=[G,g];
    end
    
    param4.G=G;
    param4.verbose=param.verbose-1;
    sol=prox_sumg(x,gamma,param4);
    
end

%reshape the solution
sol=reshape(sol,t);

iter=0;
crit='--';
info.algo=mfilename;
info.iter=iter;
info.final_eval=norm_linf1(sol,param.g_d,param.g_t);
info.crit=crit;
info.time=toc(t1);

end
