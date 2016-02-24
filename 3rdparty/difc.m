function Y = difc(varargin)
    %% circular diff function
    % difc(X,{N,D}) is a wrapper for matlab's diff()
    %
    % returns same matrix size, unlike diff().
    %
    % see diff, gradient
    %
    % created by Peter Adany
    
    %% default N,D values:
    N = 1;
    if( isvector(varargin{1}) )
        D = 1 + double( size(varargin{1},2) > size(varargin{1},1) );    % input is a vector
    else
        if( ndims(varargin{1})==2 )
            D = 1;                                                      % input is a matrix
        else
            D = min(find( size(varargin{1}) > 1 ));                     % input is an N-dimensional array
        end
    end
    
    %% accept specific N,D values:
    switch nargin
        case 1
            [];
        case 2
            N = varargin{2};
        case 3
            N = varargin{2};
            D = varargin{3};
    end
    
    %% get result using diff:
    switch N
        case 0
            Y = varargin{1};            % 0th order result is only a feature of difc, not diff
        otherwise
            Y = diff(varargin{1},1,D);  % call matlab's diff
            Y = cat(D,Y,-sum(Y,D));     % the end-wrapped first-order difference is equal to -sum(diff)
            if(N>1)
                Y = difc(Y,N-1,D);      % recursion takes care of higher order N (note, negative N is blocked by diff).
            end
    end
end