function s = generalized_forward_backward_alg()
   s.name = 'GENERALIZED_FORWARD_BACKWARD';
   s.initialize = @(x_0, fg, Fp, param) generalized_forward_backward_initialize(x_0,fg,Fp,param);
   s.algorithm = @(x_0, fg, Fp, sol, s, param) generalized_forward_backward_algorithm(fg, Fp, sol, s, param);
   s.finalize = @(x_0, fg, Fp, sol, s, param) sol;

end

function [sol, s, param] = generalized_forward_backward_initialize(x_0,fg,Fp,param)

	if ~isfield(param, 'weights'), param.weights=ones(length(Fp),1); end
    if ~isfield(param, 'lambda'), param.lambda=1 ; end
    
    % Normalizing the weights
    s.weights= param.weights./sum(param.weights);
    s.lambda = param.lambda;
    
    s.x_n = cell(length(Fp),1);
    s.u_n = cell(length(Fp),1);
    for ii = 1 : length(Fp)
        s.u_n{ii}=x_0;
    end

    sol = x_0;
    if numel(Fp)<1
        error('This solver need non smooth functions')
    end
    
    if ~fg.beta
        error('Beta = 0! This solver requires a smooth term.');
    end
    
end


function [sol, s] = generalized_forward_backward_algorithm(fg, Fp, sol, s, param)

    temp_grad = fg.grad(sol);
	for ii = 1:length(Fp)
%         s.x_n{ii} = Fp{ii}.prox_ad( 2*sol - s.u_n{ii} ...
%             - s.lambda * temp_grad ,...
%             1/s.weights(ii) * s.lambda);
%         s.u_n{ii} = s.u_n{ii} + param.gamma*(s.x_n{ii}{1} -sol);
        s.x_n{ii} = Fp{ii}.prox( 2*sol - s.u_n{ii} ...
            - s.lambda * temp_grad ,...
            1/s.weights(ii) * s.lambda);
        s.u_n{ii} = s.u_n{ii} + param.gamma*(s.x_n{ii} -sol);
	end
    
    sol = 0;
    for ii = 1:length(Fp)
        sol=sol + s.weights(ii) * s.u_n{ii};
    end

end
