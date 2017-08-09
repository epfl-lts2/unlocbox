function [im]=cameraman()
%CAMERAMAN  Load the 'cameraman' test signal
%
%   `cameraman` loads the 'cameraman' signal. The Cameraman (a.k.a.
%   Photographer) is an image commonly used in image processing, especially
%   filtering papers.  The resolution is (256 x 256).
%
%   This signal, and other standard image tests signals, can be found on
%   Morgan McGuire's Computer Graphics
%   Archive`<http://graphics.cs.williams.edu/data/images.xml>`_. 
%
%   For convenience the output image is normalized by 255 and converted to
%   double.
%
%   Example
%   -------
%   
%   Load the image and display it:::
%
%       im = cameraman();
%       imagescgray(im);
%

% Author: Nathanael Perraudin
% Date: 25 November 2013
  
if nargin>0
  error('This function does not take input arguments.')
end;

f=mfilename('fullpath');

% Load the signal
im = double(imread([f,'.png']))/255;

