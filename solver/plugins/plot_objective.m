function [ sol ] = plot_objective(x, fig)
%PLOT_OBJECTIVE Plot objective function over iters(plugin for the UNLocBoX)
%   Usage [ im ] = plot_image( im, fig );
%
%   Input parameters:
%         info_iter   : Structure of info of current iter of algorithm
%         fig   : Figure
%         setlim : Set limite for the image (default 0)
%
%   Output parameters:
%         im    : Input image
%
%   This plugin display the image every iterations of an algorithm. To use
%   the plugin juste define::
%       
%       fig=figure(100);
%       param.do_sol=@(x) plot_image(x,fig);
%
%   In the structure of optional argument of the solver.
%

% Author: Vassilis Kalofolias
% Date  : November 2013


% select the figure
figure(fig);

% display the image
subplot 211
imshow(x.sol, []);
title(['Current it: ', num2str(x.iter),'   Curr obj: ', ...
    num2str(x.curr_norm)]);
subplot 212
semilogy(x.objective); title('Objective function')
drawnow;

sol = x.sol;

end

