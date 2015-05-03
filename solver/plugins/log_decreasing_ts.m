function gamma = log_decreasing_ts(x, gamma_in, gamma_fin, nit)
%LOG_DECREASING_TS Log decreasing timestep for UNLCOBOX algorithm
%   Usage gamma = log_decreasing_ts(x, gamma_in, gamma_fin, nit);
%
%   Input parameters:
%         x         : Structure of data
%         gamma_in  : Initial timestep
%         gamma_fin : Final timestep
%         nit       : Number of iteration for the decrease
%
%   Output parameters:
%         gamma     : Timestep at iteration t
%
%   This plug-in computes a new timestep at each iteration. It makes a log
%   decreasing timestep from *gamma_in* to *gamma_fin* in *nit* iterations.
%   To use this plugin, define::
%
%       param.do_ts = @(x) log_decreasing_ts(x, gamma_in, gamma_fin, nit);
%
%   in the structure of optional argument of the solver.
%

% Author: Nathanael Perraudin
% Date:   3 April 2014


    if x.iter > nit
        gamma = gamma_fin;
    else
        ts = gamma_in./linspace(1,gamma_in/gamma_fin,nit);
        gamma = ts(x.iter);
    end
end
