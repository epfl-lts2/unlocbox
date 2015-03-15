function x = prox_tv1d(y, gamma)

% EN
% This function computes proximal operator of 1D tv norm, i.e. it solves this problem:
%    prox(x) := min_y 1/2 ||x - y||_2^2 + gamma ||y||_TV.
% Input parameters are:
%     y ... input parameter
%     gamma ... constant at tv norm
%  Use:  x = prox_tv1d(y, gamma)
%
% podle èlánku / according to article: CONDAT, Laurent. A Direct Algorithm for 1-D Total Variation Denoising. 
% IEEE Signal Processing Letters. 2013, vol. 20, issue 11. DOI: 10.1109/LSP.2013.2278339.
    
N = length(y);
x = zeros(size(y));

% a)
k = 1; k0 = 1; kplus = 1; kminus = 1; 
vmin = y(1) - gamma; vmax = y(1) + gamma;
umin = gamma; umax = - gamma;

while 1
    % b)
    if k == N
        x(N) = vmin + umin;
        break;
    end
    
    % b1)
    if y(k+1) + umin < vmin - gamma
        for i = k0:kminus
            x(i) = vmin;
        end
        k = kminus + 1; k0 = kminus + 1; kplus = kminus + 1; kminus = kminus + 1;
        vmin = y(k);
        vmax = y(k) + 2*gamma;
        umin = gamma;
        umax = - gamma;
    else
        % b2)
        if y(k+1) + umax > vmax + gamma
            for i = k0:kplus
                x(i) = vmax;
            end
            k = kplus + 1; k0 = kplus + 1; kminus = kplus + 1; kplus = kplus + 1;
            vmin = y(k) - 2*gamma;
            vmax = y(k);
            umin = gamma;
            umax = - gamma;
        else
            % b3)
            k = k + 1;
            umin = umin + y(k) - vmin;
            umax = umax + y(k) - vmax;
            % b31)
            if umin >= gamma
                vmin = vmin + (umin - gamma)/(k - k0 + 1);
                umin = gamma;
                kminus = k;
            end
            % b32)
            if umax <= -gamma
                vmax = vmax + (umax + gamma)/(k - k0 + 1);
                umax = -gamma;
                kplus = k;
            end
        end
    end
    % c)
    if k >= N
        % c1)
        if umin < 0
            for i = k0:kminus
                x(i) = vmin;
            end
            k = kminus + 1; k0 = kminus +1; kminus = kminus + 1;
            vmin = y(k);
            umin = gamma;
            umax = y(k) + gamma - vmax;
        else
            % c2)
            if umax > 0
                for i = k0:kminus
                    x(i) = vmax;
                end
                k = kplus + 1;   k0 = kplus + 1; kplus = kplus + 1;
                vmax = y(k);
                umax = - gamma;
                umin = y(k) - gamma - vmin;
            else
                % c3)
                for i = k0:N
                    x(i) = vmin + umin/(k - k0 + 1);
                end
                break;
            end
        end
    end
end

end
