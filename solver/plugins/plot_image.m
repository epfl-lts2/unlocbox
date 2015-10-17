function [ im ] = plot_image( x,fig,setlim )
%PLOT_IMAGE Plot image plugin for the UNLocBoX
%   Usage [ im ] = plot_image( im,fig );
%
%   Input parameters:
%         x     : Structure of data
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

% Author: Nathanael Perraudin
% Date  : 23 septembre 2013

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
imagesc(sol);
colormap gray;
title(['Current it: ', num2str(x.iter),'   Curr obj: ', ...
    num2str(x.info.objective(x.iter))]);
axis off         
axis image  
drawnow;
% return the image
if setlim
   x.sol = sol; 
end
im=x.sol;


end

