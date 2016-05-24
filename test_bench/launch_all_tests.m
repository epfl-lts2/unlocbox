%
%   In this function all tests are launched
%   Perraudin Nathanael, 27 nov 2012

%% Initialization

clear all;
close all;
clc;

addpath ~/work/Matlab/
%init_path

ltfatstart
init_unlocbox
%sgwt_setpath

demos=0; % to test also the demos

%% Lauch different test

errors=0;
% This need to be redone
test_proj_box; 

% Mixed sparsity
test_proj_box; 

% Prox L2
errors=errors+test_prox_l2();

% Test ppxa and other
errors=errors+test_ppxa_and_other();
errors=errors+test_ppxa();

% Test prox_lp
errors=errors+test_prox_l2grad();

% Test proj_b2
errors=errors+test_proj_b2();

% Test proj_b1
errors=errors+test_proj_b1();



% Test prox l1
errors=errors+test_prox_l1();

% Test TV
errors=errors+test_tv();

% Test solvers
errors=errors+test_solvers();

% Test prox_lp
errors=errors+test_proxlp();

% Test proj nuclear norm
errors=errors+test_proj_nuclear();

errors = errors + test_prox_functions();
errors = errors + test_utils();
errors = errors + test_lp();
errors = errors + test_prox_sum_log();
errors = errors + test_proj_simplex();
errors = errors + test_squareform();


%% Error in testbench

if errors
    fprintf('\n ****************************************** \n   They are %i test(s) ended with errors!\n ****************************************** \n',errors);
else
    fprintf('\n ****************************************** \n   No errors in testbenches!\n ****************************************** \n');
end

%% Lanch the demos
if demos
    demo_unlocbox
    demo_admm
    demo_compress_sensing
    demo_compress_sensing2
    demo_compress_sensing3
    demo_compress_sensing4
    demo_deconvolution
    demo_douglas_rachford
    demo_forward_backward
    demo_graph_reconstruction
    demo_overlaping_groups_structured_sparsity

    demo_sdmm
    demo_sound_reconstruction
    demo_tvdn
    demo_weighted_l1
end
