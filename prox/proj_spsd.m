function [sol, info] = proj_spsd( x,~,param )
%PROJ_SPSD Projection on the Symetric positive semi definite set of matrices
%   Usage:  sol=proj_spsd(x)
%           sol=proj_spsd(x, 0, param)
%           [sol,info]=proj_spsd(...)
%
%   Input parameters:
%         x     : Input signal.
%         param : Structure of optional parameters.
%   Output parameters:
%         sol   : Solution.
%         info  : Structure summarizing informations at convergence
%
%   `proj_spsd(x, gamma, param)` solves:
%
%   .. sol = argmin_{z} 0.5*||x - z||_2^2 s.t. z is SPSD
%
%   .. math::  sol = arg\min_{z} \frac{1}{2} \|x - z\|_2^2 \text{ s. t. } x \text{ is SDSD}
%   
%   param is a Matlab structure containing the following fields:
%
%   * *param.verbose* : 0 no log, 1 a summary at convergence, 2 print main
%     steps (default: 1)
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
%   * *info.crit* : Stopping criterion used 
%
%
%   See also:  prox_nuclearnorm


% Authors: Nathanael Perraudin
% Date: December 2014
%

% Start the time counter
t1 = tic;

if nargin < 3, param = struct; end

% Optional input arguments
if ~isfield(param, 'verbose'), param.verbose = 1; end



% 1) make it symetric
sol = ( x +x')/2;

% 2) make semi positive
[V,D] = eig(sol);
D(D<0) = 0;
sol = V*D*V';

% set the information structure returned
iter = 0;
crit = '--';
info.algo = mfilename;
info.iter = iter;
info.final_eval = 0;
info.crit = crit;
info.time = toc(t1);

end


