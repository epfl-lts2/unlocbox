function curr_norm = eval_function(fg,Fp,x,s,param)
%EVAL_FUNCTION internal evaluation function

% 
% if nargin<5
%     param = struct;
%     param.fast_eval = 0;
% end

curr_norm = fg.eval(x); 


for ii = 1:length(Fp)
    
%     % Here we try to accelerate the evaluation of the function by using the
%     % results of the proximal operator. The proximal functions of the
%     % UNLocBoX returns two arguments. In the second one, final_eval can be
%     % used as an approximation of the norm.
%     if param.fast_eval ...
%       && isfield(s,'x_n') ...
%       && iscell(s.x_n) ...
%       && length(s.x_n) >= ii ...
%       && iscell(s.x_n{ii}) ...
%       && length(s.x_n{ii})>1 ...
%       && isfield(s.x_n{ii}{2},'final_eval')
%         curr_norm = curr_norm + s.x_n{ii}{2}.final_eval/param.gamma;
%     else 
        curr_norm = curr_norm + Fp{ii}.eval(x); 
%     end
end

end