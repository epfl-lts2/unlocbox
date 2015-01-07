%% TEST stopstruct

tic ;         % We will measure elapsed time in a loop
             % Set up the stop box:
global FS;
        
FS = stopstruct('To stop the time','clic here') ;
% Display elapsed time
fprintf('\nSTOPSTRUCT: elapsed iterations: %5.2f\n',0)
% start the loop
for ii = 1:10000     % Check if the loop has to be stopped
   fprintf('%c',repmat(8,6,1)) ;   % clear up previous time
   fprintf('%5.0f\n',ii) ;        % display elapsed time
   
   if FS.stop()
        break;
   end
end
clear FS ;    % this structure has no use anymore
