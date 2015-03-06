function bool = kbstop(cmd)
%KBSTOP Helping function to stop loop with ctrl+d
%   Usage: kbstop('init') 
%          bool = kbstop  
%          kbstop('stop')
%          kbstop('lauched')
%
%   kbstop('init')        : initialise listening
%   bool = kbstop         : Has the caracter ctrl+d be taped ?
%   kbstop('stop')        : turn off listening
%   kbstop('lauched')     : is kbstop lauched
%
%   NOTE: There's a chance that pressing ctrl+c may interrupt the callback
%   function that saves the key being pressed. You can restart the function
%   using kbhit('init')
%
% 	Example:
%
%         tic ;         % We will measure elapsed time in a loop        
%         kbstop('init');
%         % Display elapsed time
%         fprintf('\nSTOPSTRUCT: elapsed time: %5.2f\n',toc)
%         % start the loop
%         while toc < 20    % Check if the loop has to be stopped
%             fprintf('%c',repmat(8,6,1)) ;   % clear up previous time
%             fprintf('%5.2f\n',toc) ;        % display elapsed iteration
%             if kbstop()
%                 break;
%             end
%         end
%         kbstop('stop');
%

% Author: Nathanael Perraudin
% Date  : 17 sept 2014
% Inspired by Amanda Ng of the function kbhit on matlabcentral
%
% Big acknowledgements to Yair Altman in the MATLAB newsreader thread
% "java: add keyListener to command window", for providing the exact code
% needed to set up the key press callback.  

    
global KBSTOP_h_cw
global KBSTOP_h_stop
global KBSTOP_h_cw_cbp
global KBSTOP_test_event

        
bool = 0;
drawnow;


if nargin == 0
    cmd = 'test';
end

switch lower(cmd)
    case 'test'
         if isempty(KBSTOP_h_stop)
            KBSTOP_h_stop = 0;
            bool = KBSTOP_h_stop;
         else
            bool = KBSTOP_h_stop;
         end
    case 'init'
        try %#ok<TRYNC>
            KBSTOP_test_event = @(x,y) test_event(x,y);
            KBSTOP_h_stop = 0;
            mde = com.mathworks.mde.desk.MLDesktop.getInstance;
            KBSTOP_h_cw = mde.getClient('Command Window');
            xCmdWndView = KBSTOP_h_cw.getComponent(0).getViewport.getComponent(0);
            KBSTOP_h_cw_cbp = handle(xCmdWndView,'CallbackProperties');
            javastr = ['global KBSTOP_h_stop; ' ...
                       'global KBSTOP_test_event; ' ...
                       'event = get(gcbo, ''KeyPressedCallbackData''); '...
                       'KBSTOP_h_stop = KBSTOP_test_event(event,KBSTOP_h_stop);',...
                       'clear KBSTOP_h_stop;' ... 
                       'clear KBSTOP_test_event;' ... 
                       ];
            set(KBSTOP_h_cw_cbp, 'KeyPressedCallback', javastr);
            KBSTOP_h_stop = 0;
        end
    case 'stop'
        try %#ok<TRYNC>
            set(KBSTOP_h_cw_cbp, 'KeyPressedCallback', []);
        end
        clear global KBSTOP_h_cw KBSTOP_h_stop KBSTOP_h_cw_cbp KBSTOP_test_event
    case 'lauched'
        bool = ~isempty(KBSTOP_h_stop) && ~isempty(KBSTOP_h_cw) && ...
               ~isempty(KBSTOP_h_cw_cbp) && ~isempty(KBSTOP_test_event);
        
    otherwise
        error('Unrecognised parameter');
end

end


function bool = test_event(event,bool_old)
    modifier = get(event,'modifiers');
    keycode = get(event,'KeyCode');
    bool = (modifier == 2) && (keycode == 68);
    bool = bool_old || bool;
end
