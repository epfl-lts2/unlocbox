function close_unlocbox()
%CLOSE_UNLOCBOX Closes the toolbox
%   Usage: close_unlocbox()
%
%   Close script to stop the unlocbox, release memory if gpu was used


% Author: Nathanael Perraudin
% Date: nov 2012


%% adding dependency
global GLOBAL_useGPU;

if GLOBAL_useGPU 
    reset(gpuDevice(1));
end


kbstop('stop');
