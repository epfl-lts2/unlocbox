function [im] = checkerboard()
%CHECKERBOARD Load the 'checkerboard' test signal
%   Usage: im = checkerboard();
%
%   Input parameters:
%       none
%   Output parameters:
%       im    : image
%
%   Example
%   -------
%   
%   Load the image and display it:::
%
%       im = checkerboard();
%       imagescgray(im);
%

% Author: Nathanael Perraudin
% Date: 6 February 2015
 

f = mfilename('fullpath');

% Load the signal

im = imread([f, '.png']);

im = rgb2gray(im);

im = double(im) / 255;
