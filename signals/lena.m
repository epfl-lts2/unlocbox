function [im] = lena(color)
%LENA  Load the 'lena' test signal
%   Usage: im = lena();
%          im = lena(color);
%
%   Input parameters:
%       color   : boolean 
%   Output parameters:
%       none    :
%
%   `lena()` loads the graylevel 'lena' signal. Lena is a common
%   image processing test image of resolution (512 x 512). However, we do
%   recommand not to use it. 
%
%   `lena(1)` loads the color 'lena' signal.
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
%       im = lena();
%       imagescgray(im);
%

% Author: Nathanael Perraudin
% Date: 25 November 2013
  
if nargin < 1
  color = 0;
end;

f = mfilename('fullpath');

% Load the signal

im = imread([f, '.png']);

if ~ color
    im = rgb2gray(im);
end

im = double(im) / 255;
