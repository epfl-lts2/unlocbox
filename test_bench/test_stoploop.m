%% TEST stoploop Ctrl+d


tic ;         % We will measure elapsed time in a loop
             % Set up the stop box:
        
kbstop('init');
% Display elapsed time
fprintf('\nSTOPSTRUCT: elapsed time: %5.2f\n',toc)
% start the loop
while toc < 20    % Check if the loop has to be stopped
    fprintf('%c',repmat(8,6,1)) ;   % clear up previous time
    fprintf('%5.2f\n',toc) ;        % display elapsed iteration
    if kbstop()
        break;
    end
end
kbstop('stop');
