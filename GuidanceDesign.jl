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

function circular_d(Rc, R)
    d = R - Rc;
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
        θ = atan(y/x);
        R = norm([x,y]);
        if path == "straight_line"
            d = straight_line_d(R, θ, psi_d);
        elseif path == "circular"
            d = circular_d(Rc, R)  # Example radius for circular path
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



