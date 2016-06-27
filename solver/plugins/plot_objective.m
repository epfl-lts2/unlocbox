function [ sol ] = plot_objective(info_iter, fig)
%PLOT_OBJECTIVE Plot objective function over iters(plugin for the UNLocBoX)
%   Usage [ sol ] = plot_objective( info_iter, fig );
%
%   Input parameters:
%         info_iter   : Structure of info of current iter of algorithm
%         fig   : Figure
%
%   Output parameters:
%         sol   : Current solution
%
%   This plugin displays the image every iterations of an algorithm. To use
%   the plugin juste define::
%       
%       fig = figure(100);
%       param.do_sol = @(x) plot_objective(x, fig);
%
%   In the structure of optional argument of the solver.
%

% Author: Vassilis Kalofolias
% Date  : November 2013




% select the figure
if info_iter.iter<2
    figure(fig);
end

% 
title(['Current it: ', num2str(info_iter.iter),'   Curr obj: ', ...
    num2str(info_iter.info.objective(end))]);
semilogy(info_iter.info.objective); title('Objective function')
drawnow;

% return the image
sol=info_iter.sol;

end

