function [ ] = setfig( param )
%SETFIG Set default parameters for plotting
%   Usage: setfig(param);
%   
%   Input parameters:
%       param   : optional parameters
%   Output parameters:
%       none
%
%   *param* a Matlab structure containing the following fields:
%
%   * *param.position* : position and size of the figure 
%     (default [100 100 600 400])
%   * *param.labelsize* : Size of the label (default 12)




% Nathanael Perraudin
% 25 November 2013

if nargin<1
    param=struct;
end

% Optional parameters
if ~isfield(param, 'position'), param.position = [100 100 250 250]; end
if ~isfield(param, 'labelsize'), param.labelsize = 12; end



% set the axes
   
%    set(0,'DefaultFigurePaperPosition',param.position)
    set(0,'DefaultFigurePosition',param.position);
    %set(0,'DefaultFigurePaperPosition','auto');
   
    % Change default axes fonts.
    set(0,'DefaultAxesFontSize',  param.labelsize)

    % Change default text fonts.
    set(0,'DefaultTextFontSize', param.labelsize)
    


end

