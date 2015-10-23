function [stop,rel_norm,prev_norm,iter,objectiv,crit] = convergence_test_old(curr_norm,prev_norm,iter,objectiv,param,s)
%CONVERGENCE_TEST Test the convergence of an algorithm
%   Usage: stop = convergence_test(curr_norm,prev_norm,objectiv,param);
%          stop = convergence_test(curr_norm,prev_norm,[],param);
%          stop = convergence_test(curr_norm,prev_norm);
%          [stop,rel_norm,prev_norm,iter,objectiv,crit] = convergence_test(...);
%          [~,rel_norm,prev_norm,iter,objectiv,crit] = convergence_test(curr_norm);
%          [~,rel_norm,prev_norm,iter,objectiv,crit] = convergence_test();
%
%   Input parameters:
%         curr_norm : Current norm (scalar)
%         prev_norm : Previous norm (scalar)
%         iter      : Current iteration (positiv integer)
%         objectiv  : Vector (previous objectiv function values)
%         param     : Optional parameter
%
%   Output parameters:
%         stop      : Convergence (boolean)
%         rel_norm  : Relativ norm (scalar)
%         prev_norm : Updated previous norm (scalar)
%         objectiv  : Vector (updated objectiv function values)
%         crit      : Convergence criterion
%
%   This function test the convergence of an algorithms and performs
%   updates.
%
%   At convergence, the flag stop is set to one.
%
%   If *curr_norm* and *prev_norm* are close enought (*param.tol*), the stopping
%   criterion *crit* will be 'REL_NORM'.
%
%   If the maximum number of iteration is obtained, the stoping criterion
%   *crit* is 'MAX_IT'.
%
%   If *curr_norm* is smaller than *param.tol* and *param.abs_tol* is 
%   activated, the stopping criterion *crit* will be 'CURR_NORM'.
%
%   `[~,~,prev_norm,iter,objectiv,~]= convergence_test(curr_norm)` will
%   initiate all the values for algorithm. `Prev_norm=curr_norm`
%
%   `[~,~,prev_norm,iter,objectiv,~] = convergence_test()` will
%   initiate all the values for algorithm.  `Prev_norm=eps`
%
%   *param* a Matlab structure containing the following fields:
%
%   * *param.tol* : is stop criterion for the loop. The algorithm stops if
%
%     ..  (  n(t) - n(t-1) )  / n(t) < tol,
%      
%     .. math:: \frac{  n(t) - n(t-1) }{ n(t)} < tol,
%
%     where  $n(t) = f(x)$ is the objective function at iteration *t*
%     by default, `tol=10e-4`.
%
%   * *param.maxit* : is the maximum number of iteration. By default, it is 200.
%
%   * *param.abs_tol* : If activated, this stopping critterion is the
%     objectiv function smaller than *param.tol*. By default 0.
%
%   * *param.verbose* : 0 no warning, 1 warning activated. (Default = 1)
%
%   * *param.alg* : algorithm name
%
%   * *param.alg* : Show a box to stop the algorithm (default 0)
%

%   Nathanael Perraudin
%   Date: 14 dec 2012

global FS

% Optional input arguments

% param
if nargin<5, param=struct; end
if ~isfield(param, 'tol'), param.tol=10e-4 ; end
if ~isfield(param, 'maxit'), param.maxit=200; end
if ~isfield(param, 'abs_tol'), param.abs_tol=0 ; end
if ~isfield(param, 'verbose'), param.verbose=1 ; end
if ~isfield(param, 'use_dual'), param.use_dual=0 ; end

if ~isfield(param, 'alg'), 
    txt = dbstack(); 
    param.alg = txt(2).name ; 
end

if ~isfield(param, 'stop_box'), param.stop_box = 0; end

% objectiv
if nargin<4, objectiv=[]; end

% objectiv
if nargin<3, iter=0; end

% prev_norm
if nargin<2, prev_norm=0; end

% curr_norm
if nargin<1
   curr_norm=0;
   param.verbose=0;
end

if ~kbstop('lauched')
    kbstop('init');
end

if param.use_dual &&  nargin>=6 && isfield(s,'reldual')
    curr_norm = s.reldual;
    param.abs_tol = 1;
end
if param.stop_box
    if iter<=1
        if isstruct(FS)
            try  %#ok<TRYNC>
                FS.close();
            end
        end
        FS = stopstruct(param.alg,'Stop the algorithm') ;
    end
end
% % perform simple test
% if (curr_norm==0)
%     if param.verbose
%         fprintf('WARNING: current norm is equal to 0! Adding eps to continure...\n');
%     end
%     curr_norm=eps;
% end

if ~(numel(curr_norm)==1)
    error('One of your evaluation function does not return a scalar.')
end

rel_norm = abs(curr_norm - prev_norm)/(curr_norm + eps);
if iter
    if isa(curr_norm,'gpuArray')
        objectiv(iter)=gather(curr_norm);
    else
        objectiv(iter)=curr_norm;
    end
end

if numel(curr_norm)>1
    error('One of your evaluation functions .eval returns a vector and not a scalar');
end

if (curr_norm < param.tol) && logical(param.abs_tol) && iter>1
    crit = 'CURR_NORM';
    stop=1;
elseif (abs(rel_norm) < param.tol) && (~param.abs_tol) && iter>1
    crit = 'REL_NORM';
    stop=1;
elseif iter >= param.maxit
    crit= 'MAX_IT';
    stop=1;    
% would be great to have a text like this (thinking about it)    
% elseif curr_norm > prev_norm
%     crit = 'OBJ_INC';
%     stop = 1;
elseif  kbstop() || (param.stop_box && FS.stop())
    crit= 'USER';
    stop = 1;
else
    crit= 'NOT_DEFINED';
    stop=0;
end

% Performed updates
if ~stop
    iter=iter+1;
end
prev_norm=curr_norm;

if stop
    kbstop('stop');
    if param.stop_box

        try %#ok<TRYNC>
            FS.close() ;  % Clear up the box
        end

        try %#ok<TRYNC>
            clear FS ;    % this structure has no use anymore 
        end
    end
end