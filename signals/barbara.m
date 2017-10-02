function [im]=barbara()
%BARBARA  Load the 'barbara' test signal
%
%   `barbara` loads the 'barbara' signal. Barbara is an image commonly
%   used in image compression and filtering papers because it contains a
%   range of tones and many thin line patterns. The resolution is (512 x
%   512). 
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
%       im = barbara();
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

