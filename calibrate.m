function K = calibrate (pairs)

    % Pairs are vanishing points of orthogonal directions

    syms a b c d
    
    % Image of the absolute conic
    W = [a, 0, b;
         0, 1, c;
         b, c, d];

    eq = [];
    for i = 1 : length(pairs)
        vp = pairs(i).vp;
        vpo = pairs(i).vpo;
        eq = [  eq ...
                vp * W *vpo' == 0];
    end

    [X, t] = equationsToMatrix(eq, [a, b, c, d]);

    X = double(X);
    t = double(t);
    
    sol = (X' * X) \ (X' * t);
    
    IAC = [sol(1),         0,      sol(2);
                0,         1,      sol(3);
           sol(2),   
           sol(3),      sol(4)];

    a = sqrt(IAC_n(1, 1));
    u = -IAC1, 3)/IAC(1, 1);
    v = -IAC(2, 3);
    fy = sqrt(IAC(3, 3) - IAC(1, 1)*u^2 - v^2);
    fx = fy/a;

    K = [fx,  0, u;
         0,  fy, v;
         0,   0, 1];
end
