function [ s ] = plot_signal( x,p,fig )
%PLOT_SIGNAL Plot image plugin for the UNLocBoX
%   Usage [ s ] = plot_signal( im,fig );
%
%   Input parameters:
%         x     : Structure of data
%         p     : Number of iteration between 2 plots...
%         fig   : Figure
%
%   Output parameters:
%         s     : Input image
%
%   This plugin display the signal every iterations of an algorithm. To use
%   the plugin juste define::
%       
%       fig=figure(100);
%       param.do_sol=@(x) plot_signal(x,p,fig);
%
%   In the structure of optional argument of the solver.
%

% Author: Nathanael Perraudin
% Date  : 3rd april 2014

if ~mod(x.iter-1,p)
    % select the figure
    if x.iter<2
        figure(fig);
    end
    % display the signal
    plot(x.sol);
    title(['Current it: ', num2str(x.iter),'   Curr obj: ', ...
        num2str(x.info.objective(x.iter))]);

    drawnow;
end

% return the signal
s=x.sol;


end

