% Unlocbox - Proximal operators
%
%   The proximal operator of function f evaluated in z is the solution of 
%   the folowing equation:   
%
%   .. prox_{f, gamma }(x)=min_x  1/2 ||x-z||_2^2 + gamma f(z)
%
%   .. math:: prox_{f, \gamma }(z)=\min_z  \frac{1}{2} \|x-z\|_2^2 + \gamma f(z)
%
%   Here are a list of common usual proximal operators available in the
%   UnLocBoX. We remember the reader that projections are particular cases of proximal
%   operators.
%
%  General Proximal operators
%    PROX_L0            -  Proximal operator of the L0 norm
%    PROX_L1            -  Proximal operator of the L1 norm
%    PROX_L2            -  Proximal operator of the L2 norm
%    PROX_L2grad        -  Proximal operator of the L2 norm of the gradient 
%    PROX_L2gradfourier -  Proximal operator of the L2 norm of the gradient in the Fourier domain
%    PROX_Linf1         -  Proximal operator of the Linf1 norm
%    PROX_L21           -  Proximal operator of the L21 norm
%    PROX_L12           -  Proximal operator of the LP norm
%    PROX_NUCLEARNORM   -  Proximal operator of the nuclear norm
%    PROX_NUCLEARNORM_BLOCK -  Proximal operator of the nuclear norm by block
%    PROX_TV            -  Proximal operator of the TV norm
%    PROX_TV3D          -  Proximal operator of the 3D TV norm
%    PROX_TV1D          -  Proximal operator of the 1D TV norm
%    PROX_TV4D          -  Proximal operator of the 4D TV norm
%    PROX_SUM_LOG       -  Proximal operor for the sum of logarithms.
%    PROX_SUM_LOG_NORM2 -  Proximal operator of log-barrier  - sum(log(x))
%
%  Projection operators
%    PROJ_B1            -  Projection on a B1-Ball
%    PROJ_B2            -  Projection on a B2-Ball
%    PROJ_BOX           -  Bound the value of x
%    PROJ_NUCLEARNORM   -  Projection operator of the nuclear norm
%    PROJ_SPSD          -  Projection onto the semi definite symetric positive set 
%    PROJ_LINEAR_EQ     -  Projection onto a linear set of equalities
%    PROJ_LINEAR_INEQ   -  Projection onto a linear set of inequalities
%   
%
%  Proximal tools
%    PROX_SUMG          -  Proximal operator of a sum of function
%    PROX_ADJOINT       -  Proximal operator of the adjoint of a function
%    PROX_ADD_2NORM     -  Proximal operator
%    PROX_FAX           -  Proximal operator of f(Ax)
%
%  For help, bug reports, suggestions etc. please send email to
%  unlocbox (at) groupes (dot) epfl (dot) ch
%
% see also: solver solvep

