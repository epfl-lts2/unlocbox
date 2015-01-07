function [ errors ] = test_ppxa( )
%TEST_PROX_L2 This function test the function prox_l2
errors=0;

errors=errors+test_ppxa_only(10e-5);

end

function [errors]=test_ppxa_only(tol)

     [x0, f1,f2,param ] = create_cs_problem();
     param.verbose=0;
     param.maxit=5;
     p2=ppxa_old(x0,[f1,f2],param);
     p3=ppxa(x0,{f1,f2},param);
    if norm(p3-p2)/norm(p3)<tol
        fprintf('  Test ppxa OK\n')
        errors=0;
    else
        fprintf('  Test ppxa Pas OK!!!!!!!!!!!!!!!!\n')
        norm(p3-p2)/norm(p3)
        errors=1;
    end
end
