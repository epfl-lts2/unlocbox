function [sol,param]=post_process(sol,iter,curr_norm,prev_norm,param)
%POST_PROCESS Post processing for the UNLocBoX
%   Usage: [sol,param]=post_process(sol,iter,curr_norm,prev_norm,param);
%
%   
%   This function make evaluate the post processing job in the UnLocBoX

if ~isfield(param, 'gamma'), param.gamma=1; end;
if ~isfield(param, 'do_sol'), param.do_sol=@(x) x.sol ; end
if ~isfield(param, 'do_ts'), param.do_ts=@(x) x.gamma ; end

x.sol=sol;
x.iter=iter;
x.curr_norm=curr_norm;
x.prev_norm=prev_norm;
x.gamma=param.gamma;

sol=param.do_sol(x);
param.gamma=param.do_ts(x);



end

