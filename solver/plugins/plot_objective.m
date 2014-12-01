function [ sol ] = plot_objective(x, fig,setlim)
%PLOT_OBJECTIVE Plot objective function over iters(plugin for the UNLocBoX)
%   Usage [ im ] = plot_image( im, fig );
%         [ im ] = plot_image( im, fig , setlim);
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


if nargin<3
    setlim = 0; 
end

% select the figure
if x.iter<2
    figure(fig);
end

sol = x.sol;
sol(sol<0) = 0;
sol(sol>1) = 1;

% display the image
subplot 211
imshow(sol, []);
title(['Current it: ', num2str(x.iter),'   Curr obj: ', ...
    num2str(x.curr_norm)]);
subplot 212
semilogy(x.objective); title('Objective function')
drawnow;

% return the image
if setlim
   x.sol = sol; 
end
sol=x.sol;

end

