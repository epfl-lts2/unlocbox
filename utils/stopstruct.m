function FS = stopstruct(titletxt,buttontxt)
%STOPSTRUCT Create a stoping structure for a loop
%   Usage:  FS = stopstruct(titletxt,buttontxt);
%   
%   This function create a structure to interrupt for or while loop.
%
%   Example::
%   
%         global FS;
% 
%         FS = stopstruct('To stop the time','clic here') ;
%         % Display elapsed time
%         fprintf('\nSTOPSTRUCT: elapsed iterations: %5.2f\n',0)
%         % start the loop
%         for ii = 1:10000     % Check if the loop has to be stopped
%            fprintf('%c',repmat(8,6,1)) ;   % clear up previous time
%            fprintf('%5.0f\n',ii) ;        % display elapsed time
% 
%            if FS.stop()
%                 break;
%            end
%         end
%         clear FS ;    % this structure has no use anymore

FS.state = 0;


% make a uicontrol that will be used to interrupt some loop
% intu = 65756;
FS.fig = figure(...
    'Units','characters',...
'PaperUnits',get(0,'defaultfigurePaperUnits'),...
'IntegerHandle','off',...
'InvertHardcopy',get(0,'defaultfigureInvertHardcopy'),...
'MenuBar','none',...
'Name',titletxt,...
'NumberTitle','off',...
'PaperPosition',get(0,'defaultfigurePaperPosition'),...
'PaperSize',get(0,'defaultfigurePaperSize'),...
'PaperType',get(0,'defaultfigurePaperType'),...
'Position',[100 50 50 10],...
'Resize','off',...
'HandleVisibility','callback',...
'Visible','on');


FS.close = @() closefig(FS.fig);
FS.stop = @() statefig;


set(FS.fig,'HandleVisibility','on');



FS.UI = uicontrol(...
    'style','push',...
    'callback','global FS; FS.state = 1; FS.close();',...
    'Units','characters',...
    'FontSize',20,...
    'Position',[4 2 42 6],...
    'String',buttontxt);
winontop(FS.fig);
%title(str);
set(FS.fig,'HandleVisibility','off')




end


function closefig(fig)
    try %#ok<TRYNC>
        set(fig,'HandleVisibility','on');
        close(fig);
    end
    
end

function state = statefig()
    global FS
    drawnow;
    state = FS.state;
end
