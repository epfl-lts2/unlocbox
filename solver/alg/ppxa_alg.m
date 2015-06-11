function s = ppxa_alg()
   s.name = 'PPXA';
   s.initialize = @(x_0, fg, Fp, param) ppxa_initialize(x_0,fg,Fp,param);
   s.algorithm = @(x_0, fg, Fp, sol, s, param) ppxa_algorithm(Fp, sol, s, param);
   s.finalize = @(x_0, fg, Fp, sol, s, param) sol;

end

function [sol, s, param] = ppxa_initialize(x_0,fg,Fp,param)

    if ~isfield(param, 'lambda'), param.lambda=0.99 ; end
    if ~isfield(param, 'W'), param.W = ones(length(Fp), 1) / length(Fp); end
    s.lambda = param.lambda;
    s.W = param.W;
    
    % Reshape x if vector line
    if (size(s.W, 2) > size(s.W, 1))
        s.W = s.W';
    end
    test_sum(s.W);
    
    
    
    % Create a table of cell containing data
    s.y_n = cell(length(Fp),1);
    s.x_n = cell(length(Fp),1);

    for ii = 1:length(Fp)
        s.y_n{ii} = x_0;    
%        s.x_n{ii}{1} = zeros(size(x_0));  
        s.x_n{ii} = zeros(size(x_0));  
    end
        
    sol = x_0;
    
    if fg.beta
        error('PPXA needs only function with proximal operators')
    end
    
    
end


function [sol, s] = ppxa_algorithm(Fp, sol, s, param)

    % proximal operator
    for ii = 1:length(Fp)
%        s.x_n{ii} = Fp{ii}.prox_ad(s.y_n{ii}, param.gamma);
        s.x_n{ii} = Fp{ii}.prox(s.y_n{ii}, param.gamma);
    end
    
    pn = w_sum(s.W,s.x_n);

    % update y

    for ii = 1:length(Fp)
%        s.y_n{ii}  = s.y_n{ii} + param.lambda * (2*pn - sol - s.x_n{ii}{1}); 
        s.y_n{ii}  = s.y_n{ii} + param.lambda * (2*pn - sol - s.x_n{ii}); 
    end
    
    % update sol
    sol = sol + param.lambda * (pn - sol);

       
end


function test_sum(W)
    if (sum(W) > 1+eps) || (sum(W) < 1-eps)
        fprintf('Warning : sum W is not equal to 1');
    end
end

function s = w_sum(w, data)

if iscell(data{1})
    s = zeros(size(data{1}{1}));
    for ii = 1:length(data)
        s = s + w(ii) * data{ii}{1};
    end    
else
    s = zeros(size(data{1}));
    for ii = 1:length(data)
        s = s + w(ii) * data{ii};
    end
end

end

