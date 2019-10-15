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
%   * *param.epsilon* : Radius of the L1 ball (default = 1).
%
%   * *param.weight* : contain the weights (default ones).
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
%   This code is partly borrowed from the SPGL toolbox!
%
%   See also:  proj_b2 prox_l1

%
% Author: Nathanael Perraudin
% Date: February 2015
% Testing: test_proj_b1

% Start the time counter
t1 = tic;

% Optional input arguments
if nargin<3, param=struct; end
if ~isfield(param, 'epsilon'), param.epsilon = 1; end
if ~isfield(param, 'verbose'), param.verbose = 1; end
if ~isfield(param, 'weight'), param.weight = ones(size(x)); end

if isfield(param,'w')
    error('Change in the UNLocBoX! Use weight instead of w!');
end

if isscalar(param.weight), param.weight = ones(size(x))* param.weight; end
param.weight = abs(param.weight);

% Quick return for the easy cases.
if sum(param.weight) == 0
   sol   = x;
   iter = 0;

    crit='--';
    info.algo=mfilename;
    info.iter=iter;
    info.final_eval=0;
    info.crit=crit;
    info.time=toc(t1);
   return
end

% Get sign of b and set to absolute values
is_real = isreal(x);
if is_real
    signx = sign(x);
else
    phi = angle(x);
end
x = abs(x);

idx = find(x > eps); % Get index of all non-zero entries of d
sol   = x;             % Ensure x_i = b_i for all i not in index set idx
[sol(idx),iter] = one_projector(sol(idx),param.weight(idx),param.epsilon);


% Restore signs in x
if is_real
    sol = sol.*signx;
else
    sol = sol.*exp(1i*phi);
end


% Log after the projection onto the L2-ball
if param.verbose >= 1
    fprintf('  Proj. B1: epsilon = %e, ||x||_2 = %e,\n', param.epsilon, norm(sol,1));
end

crit='--';
info.algo=mfilename;
info.iter=iter;
info.final_eval=norm(param.weight.*sol,1);
info.crit=crit;
info.time=toc(t1);

end




function [sol,iter] = one_projector(x,weight,tau)
% This code is partly borrowed from the SPGL toolbox
  % Initialization
   N = length(x);
   sol = zeros(N,1);

   % Check for quick exit.
   if (tau >= norm(weight.*x,1)), sol = x; iter = 0; return; end
   if (tau <  eps         ),        iter = 0; return; end

   % Preprocessing (b is assumed to be >= 0)
   [sw,idx] = sort(x ./ weight,'descend'); % Descending.
   x  = x(idx);
   weight  = weight(idx);

   % Optimize
   csdb = 0; csd2 = 0;
   soft = 0; ii = 1;
   while (ii <= N)
      csdb = csdb + weight(ii).*x(ii);
      csd2 = csd2 + weight(ii).*weight(ii);
  
      alpha1 = (csdb - tau) / csd2;
      alpha2 = sw(ii);

      if alpha1 >= alpha2
         break;
      end
    
      soft = alpha1;  ii = ii + 1;
   end
   sol(idx(1:ii-1)) = x(1:ii-1) - weight(1:ii-1) * max(0,soft);

   % Set number of iterations
   iter = ii;

end
