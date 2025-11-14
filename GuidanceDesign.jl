
struct S_d
    Beta :: Float64
    p :: Int64
    q :: Int64
    eta :: Float64
end

struct GuidanceParams
    delta_max :: Float64
    S_d :: S_d
end

function compute_S_d(d, d_dot, GP :: GuidanceParams)
    S_d = d + (1/GP.S_d.Beta) * sign(d_dot)*abs(d_dot)^(GP.S_d.p / GP.S_d.q);
    return S_d
end

function compute_a_mc(psi, psi_d, d_dot, S_d_val, GP :: GuidanceParams)
    alpha = GP.S_d.p / GP.S_d.q;
    c = cos(psi - psi_d);
    # prevent blow-up if heading difference ≈ 90°
    if abs(c) < 1e-6
        c = sign(c) * 1e-6;
    end
    a_mc = - (1 / c) * (GP.S_d.Beta * alpha * abs(d_dot)^(2-alpha) + GP.S_d.eta * sign(S_d_val));
    return a_mc
end

mutable struct SimParams
    GP :: GuidanceParams
    v :: Float64
end

function straight_line_d(R, θ, psi_d)
    d = R * sin(θ - psi_d);
    return d
end

function circular_d(R, R_prime)
    d = R_prime - R;
    return d
end

function simulate!(
    SP :: SimParams, psi, psi_d, init_pos,
    path="straight_line", N=2000, dt=0.01, Rc=20)
    X = zeros(N);
    Y = zeros(N);
    x, y = init_pos;

    X[1] = x;
    Y[1] = y;
    for t = 2:N
        if path == "straight_line"
            θ = atan(y / x);
            R = norm([x,y]);
            d = straight_line_d(R, θ, psi_d);
        elseif path == "circular"
            pos = [X[t-1], Y[t-1]];
            n = normalize(pos - c);
            grad_R = (n[2]-c[2] ) / ( n[1]-c[1] );
            grad_tangent = - 1 / grad_R;
            psi_d = atan(grad_tangent, 1);
            # States

            R_prime = norm([X[t-1], Y[t-1]] - c);
            d = circular_d(Rc, R_prime);
        end
        d_dot = SP.v * sin(psi-psi_d);
        S = compute_S_d(d, d_dot, GP);
        a_mc = compute_a_mc(psi, psi_d, d_dot, S, GP);

        psi_dot = a_mc / SP.v;
        x = x + dt * SP.v * cos(psi);
        y = y + dt * SP.v * sin(psi);
        psi = psi + dt * psi_dot;
        Y[t] = y;
        X[t] = x;
    end
    return X, Y
end



