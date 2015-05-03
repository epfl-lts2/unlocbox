function init_unlocbox()
%INIT_UNLOCBOX Initialize the toolbox
%   Usage: init_unlocbox()
%
%   Initialisation script for the unlocbox
%   This script add the different path needed to run the toolbox


% Author: Nathanael Perraudin
% Date: nov 2012


%% adding dependency
global GLOBAL_useGPU;
global GLOBAL_path;
GLOBAL_path = fileparts(mfilename('fullpath'));
GLOBAL_useGPU = 0;

addpath(genpath(GLOBAL_path));

% Load the version number
bp=[GLOBAL_path,filesep];
[FID, MSG] = fopen ([bp,'unlocbox_version'],'r');
if FID == -1
    error(MSG);
else
    unlocbox_version = fgetl (FID);
    fclose(FID);
end

banner = sprintf(strcat(... 
'UnLocBoX version %s. Copyright 2012-2015 LTS2-EPFL, by Nathanael Perraudin'), ...
                   unlocbox_version);
% display banner
disp(banner);

if GLOBAL_useGPU && gpuDeviceCount
    dev=gpuDevice(1);
    if dev.DeviceSupported
        reset(gpuDevice);
        disp('GPU loaded');
    else
        disp(['GPU not loaded.  To remove the previous warning, '...
        'set GLOBAL_useGPU to 0']);  
        GLOBAL_useGPU=0;
    end
end

kbstop('stop');