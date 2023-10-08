function [y] = MM_lp(Tx, Rx, Rg, p, epsilon, Nmax_in, Nmax_out)

[H, M] = size(Tx);
[~, N] = size(Rx);

cnt = 0;

y = 200*(rand(H,1)-0.5);

while 1

    cnt = cnt + 1;

    y_old = y;

    w_mtx = zeros(M,N);

    for m = 1:M
        for n = 1:N
            w_mtx(m,n) = ((Rg(m,n) - norm(y - Tx(:,m)) - norm(y - Rx(:,n)))^2 + 0.000001)^((p-2)/2);
        end
    end
    
    
    
    %update y
    y = MM_l2(Tx, Rx, Rg, w_mtx, epsilon, Nmax_in);


    if ( (norm(y - y_old)/min(norm(y),norm(y_old))) < epsilon ) || (cnt > Nmax_out)
        break
    end
    
end


end

