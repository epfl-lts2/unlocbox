function [sol, param] = post_process(sol, iter, info, param)
%POST_PROCESS Post processing for the UNLocBoX
%   Usage: [sol,param]=post_process(sol,iter,curr_norm,prev_norm,param);
%
%   
%   This function make evaluate the post processing job in the UnLocBoX

if ~isfield(param, 'gamma'), param.gamma = 1; end;
if ~isfield(param, 'do_sol'), param.do_sol = @(x) x.sol; end
if ~isfield(param, 'do_ts'), param.do_ts = @(x) x.gamma ; end

info_iter.sol = sol;
info_iter.iter = iter;
info_iter.gamma = param.gamma;
info_iter.info = info;

% Perform tunable user actions every iterations like display the current
% solution or update the timesteps
sol = param.do_sol(info_iter);
param.gamma = param.do_ts(info_iter);



end

