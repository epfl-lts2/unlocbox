function [] = imagesc_gray( im,nfig,tit,subplotn,clim )
%IMAGESCGRAY Display an image in gray
%   Usage: imagescgray(im);
%          imagescgray(im,nfig);
%          imagescgray(im,nfig,title);
%          imagescgray(im,nfig,title,subplot);
%
%   Input parameters:
%         Im     : Image in matrix form
%         nfig   : Number of the figure 
%         tit    : Title of the image (string)
%         subplotn: Number of the subplot
%         clim   : Limit for the imagesc
%   Output parameters:
%
%   This function display an image in gray on this form::
%   
%        figure(nfig);
%        subplot(subplot)
%        imagesc(im,clim);
%        colormap gray;
%        hold on;
%        title(tit);
%        axis off;          
%        axis image;
%        drawnow;
%
%   If $nfig=0$, then function will automatically create a new figure
%   using::
%        
%        figure();
%
%   Example:::
%   
%       img = lena();
%       imagesc_gray(img);
%

% Author: Nathanael Perraudin
% Date  : 14.03.2013


if nargin > 1
    if nfig == 0
        figure();
    else
        figure(nfig);
    end
end
if nargin > 3
    subplot(subplotn)
end
if nargin > 4  
    imagesc(im, clim);
else
    imagesc(im);
end
colormap gray;
hold on;
if nargin > 2
    title(tit);
end
axis off;          
axis image;
X = get(gcf, 'PaperPosition');
set(gcf, 'PaperPosition', [ X(1) X(2) .5 * X(3) .5 * X(4) ] );

end

