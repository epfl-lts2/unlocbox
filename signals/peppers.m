function [im]=peppers(color)
%PEPPERS  Load the 'peppers' test signal
%   Usage: im = peppers();
%          im = peppers(color);
%
%   Input parameters:
%       color   : boolean 
%   Output parameters:
%       none    :
%
%   `peppers()` loads the graylevel 'peppers' signal. Peppers is a common
%   image processing test image of resolution (512 x 512).   
%
%   `peppers(1)` loads the color 'peppers' signal.
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
%       im = peppers();
%       imagescgray(im);
%  

% Author: Nathanael Perraudin
% Date: 25 November 2013
  
if nargin<1
  color = 0;
end;

f=mfilename('fullpath');

% Load the signal

im = imread([f,'.png']);

if ~color
    im = rgb2gray(im);
end

im = double(im)/255;

