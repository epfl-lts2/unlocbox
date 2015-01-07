function [sol,info] = proj_b1(x, ~, param)
%PROJ_B1 Projection onto a L1-ball
%   Usage:  sol=proj_b1(x, ~, param)
%           [sol,infos]=proj_b1(x, ~, param)
%
%   Input parameters:
%         x     : Input signal.
%         param : Structure of parameters.
%   Output parameters:
%         sol   : Solution.
%         info  : Structure summarizing informations at convergence
%
%   `proj_b1(x,~,param)` solves:
%
%   .. sol = argmin_{z} ||x - z||_2^2   s.t.  ||w.*z||_1 < epsilon
%
%   .. math::  sol = \min_z ||x - z||_2^2 \hspace{1cm} s.t. \hspace{1cm}  \|w.*z\|_1 < \epsilon
%
%   Remark: the projection is the proximal operator of the indicative function of
%   $||w.*z||_1 < \epsilon$. So it can be written:
%
%   .. prox_{f, gamma }(x)      where       f= i_c(||w.*z||_1 < epsilon)
%
%   .. math:: prox_{f, \gamma }(x) \hspace{1cm} where \hspace{1cm} f= i_c(\|w.*z\|_1 < \epsilon)
%
%   param is a Matlab structure containing the following fields:
%
%   * *param.epsilon* : Radius of the L1 ball (default = 1e-3).
%
%   * *param.w* : contain the weights (default ones).
%
%   * *param.verbose* : 0 no log, 1 a summary at convergence, 2 print main
%     steps (default: 1)
%
%
%   info is a Matlab structure containing the following fields:
%
%   * *info.algo* : Algorithm used
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
%   Rem: The input "~" is useless but needed for compatibility issue.
%
%   This function just recalls the function oneProjector from spgl toolbox.
%   It's has been done for compatibility issues.
%
%   See also:  proj_b2 prox_l1

%
% Author: Nathanael Perraudin
% Date: October 2011
%

% Start the time counter
t1 = tic;

% Optional input arguments
if nargin<3, param=struct; end
if ~isfield(param, 'epsilon'), param.epsilon = 1e-3; end
if ~isfield(param, 'verbose'), param.verbose = 1; end
if ~isfield(param, 'w'), param.w = ones(size(x)); end

% Call the function one Projector
sol= oneProjector(x,param.w,param.epsilon);

% Log after the projection onto the L2-ball
if param.verbose >= 1
    fprintf('  Proj. B1: epsilon = %e, ||x||_2 = %e,\n', param.epsilon, norm(sol,1));
end

iter=0;
crit='--';
info.algo=mfilename;
info.iter=iter;
info.final_eval=norm(sol,1);
info.crit=crit;
info.time=toc(t1);

end
