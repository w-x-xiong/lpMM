function [y] = MM_l2(Tx, Rx, Rg, w_mtx, epsilon)

[H, M] = size(Tx);
[~, N] = size(Rx);

y = 200*(rand(H,1)-0.5);

cnt = 0;

while 1

    cnt = cnt + 1;

    y_old = y;

    g_mtx = zeros(M,N);
    p_tensor = zeros(M,N,H);
    q_tensor = zeros(M,N,H);
    a_mtx = zeros(M,N);
    b_mtx = zeros(M,N);
    u_mtx = zeros(M,N);
    v_mtx = zeros(M,N);
    
    for m = 1:M
        for n = 1:N
            g_mtx(m,n) = Rg(m,n) - norm(y - Tx(:,m)) - norm(y - Rx(:,n));
            p_tensor(m,n,:) = Rg(m,n)*(y - Tx(:,m))/norm(y - Tx(:,m));
            q_tensor(m,n,:) = Rg(m,n)*(y - Rx(:,n))/norm(y - Rx(:,n));
            a_mtx(m,n) = norm(y - Rx(:,n))/norm(y - Tx(:,m));
            b_mtx(m,n) = norm(y - Tx(:,m))/norm(y - Rx(:,n));
            u_mtx(m,n) = g_mtx(m,n)/(2*norm(y - Tx(:,m)));
            v_mtx(m,n) = g_mtx(m,n)/(2*norm(y - Rx(:,n)));
        end
    end
    
    %update y
    numerator = zeros(H,1);
    denominator = 0;

    for m = 1:M
        for n = 1:N
            numerator = numerator + w_mtx(m,n)*((1+a_mtx(m,n)+u_mtx(m,n))*Tx(:,m) + (1+b_mtx(m,n)+v_mtx(m,n))*Rx(:,n) + vec(p_tensor(m,n,:)) + vec(q_tensor(m,n,:)));
            denominator = denominator + w_mtx(m,n)*(1+a_mtx(m,n)+u_mtx(m,n)+1+b_mtx(m,n)+v_mtx(m,n));
        end
    end

    y = numerator/denominator;


    if ( (norm(y - y_old)/norm(y_old)) < epsilon ) || (cnt > 1000)
        break
    end
    
end


end

