function algo = algoname(name)
%ALGONAME return an algorithm from is name
%   Usage: algo = algoname(name)
%
%   Input parameters:
%       name    : name of the algorithm (string)
%   Output parameters:
%       algo    : algorithm (struct)
%
%   The structure algo contains 3 fields:
%   * *algo.name* : the name of the algorithm (string)
%   * *algo.ignitialize* : the initialization function of the algorithm
%   * *algo.algorithm* : the core of one iteration of the algorithm
%   * *algo.finalize* : post processing
%


switch lower(name)
    case 'forward_backward'
        algo = forward_backward_alg();
    case 'douglas_rachford'
        algo = douglas_rachford_alg();
    case 'admm'
        algo = admm_alg();
    case 'sdmm'
        algo = sdmm_alg();
    case 'ppxa'
        algo = ppxa_alg();        
    case 'generalized_forward_backward'
        algo = generalized_forward_backward_alg();      
    case 'gradient_descent'
        algo = gradient_descent_alg();
%     case 'backward_backward'
%         algo = backward_backward_alg();  
    case 'pocs'
        algo = pocs_alg();
    case 'chambolle_pock'
        algo = chambolle_pock_alg();      
    case 'fb_based_primal_dual'
        algo = fb_based_primal_dual_alg();      
    case 'fbf_primal_dual'
        algo = fbf_primal_dual_alg();  
    otherwise
        error('Unknown algorithm name')
        
end

end

