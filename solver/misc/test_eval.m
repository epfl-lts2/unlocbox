function [F] = test_eval(F)
%TEST_NORM Test if the norm is given with the actual way
%   Usage: F = test_eval(F);
%          test_eval(F);
%
%   Input parameters:
%         F         : Scuctures or array of structures
%
%   Output parameters:
%         F         : Scuctures or array of structures
%
%   This function test if the user still use the old *.norm* name instead of
%   *.eval*. In this case it will fix the bug by coppying *.norm* to
%   *.eval* and display a warning.
%
%   This function also check if the function contain a *.norm* or a
%   *.eval*. If not, it returns an error.
%
%   At convergence, the flag stop is set to one.
%

%   Nathanael Perraudin
%   Date: 22 jan 2013

% number of function
m = length(F);
if m>1 || iscell(F)

    for ii=1:m
        F{ii} = test_eval(F{ii});
        if ~isa(F{ii}.eval,'function_handle')
            error('f.eval is not a funtion handle');
        end
    end
    
else
    flag=0;
    
    if ~isfield(F, 'eval')
        if ~isfield(F, 'norm')
            error('No eval function!');
        else
            flag=1;
            F.eval=F.norm;
        end
    end

    if flag
       fprintf('WARNING!!! You are using the old name "norm" to evalute a function. \n    Please upgrade to the new name: "eval". \n    This name will not be supported anymore in the future versions of the UNLocBoX.\n'); 
    end
end

