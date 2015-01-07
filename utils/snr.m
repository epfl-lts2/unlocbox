function snr_val = snr(map_init, map_noisy)
%SNR Compute the SNR between two maps
%   Usage snr_val = snr(map_init, map_noisy) 
% 
%   Input parameters:
%         map_init : initial signal
%         map_recon: noisy signal
%   Output parameters:
%         snr_val  : snr
%
%   computes the SNR between the maps `map_init`
%   and `map_noisy`. The SNR is computed as:
%
%       `10 * log10( var(map_init) / var(map_init-map_noisy) )`
%
%   where var stands for the matlab built-in function that computes the
%   variance.
% 
 
% Author: Gilles Puy
% Date: 2009
% 

noise = map_init(:)-map_noisy(:);
var_init = var(map_init(:));
var_den = var(noise(:));
snr_val = 10 * log10(var_init/var_den);

end
