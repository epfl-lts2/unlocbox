function [sol, param] = post_process_old(sol, iter, curr_eval, prev_eval, objective, param)
%POST_PROCESS Post processing for the UNLocBoX
%   Usage: [sol,param]=post_process(sol,iter,curr_norm,prev_norm,param);
%
%   
%   This function make evaluate the post processing job in the UnLocBoX

if ~isfield(param, 'gamma'), param.gamma=1; end;
if ~isfield(param, 'do_sol'), param.do_sol=@(x) x.sol ; end
if ~isfield(param, 'do_ts'), param.do_ts=@(x) x.gamma ; end

info_iter.sol = sol;
info_iter.iter = iter;
info_iter.curr_eval = curr_eval;
info_iter.prev_eval = prev_eval;
info_iter.gamma = param.gamma;
info_iter.objective = objective;

% TODO: why do we have this?
sol = param.do_sol(info_iter);
param.gamma = param.do_ts(info_iter);



end