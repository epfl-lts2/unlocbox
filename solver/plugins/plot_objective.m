function [ sol ] = plot_objective(x, fig)
%PLOT_OBJECTIVE Plot objective function over iters(plugin for the UNLocBoX)
%   Usage [ im ] = plot_image( im, fig );
%
%   Input parameters:
%         info_iter   : Structure of info of current iter of algorithm
%         fig   : Figure
%
%   Output parameters:
%         im    : Input image
%
%   This plugin display the image every iterations of an algorithm. To use
%   the plugin juste define::
%       
%       fig=figure(100);
%       param.do_sol=@(x) plot_objective(x,fig);
%
%   In the structure of optional argument of the solver.
%

% Author: Vassilis Kalofolias
% Date  : November 2013




% select the figure
if x.iter<2
    figure(fig);
end

% 
title(['Current it: ', num2str(x.iter),'   Curr obj: ', ...
    num2str(x.curr_norm)]);
semilogy(x.objective); title('Objective function')
drawnow;

% return the image
sol=x.sol;

end

