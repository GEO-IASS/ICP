function [R, T, detX, error ] = ICP_showy(P, P_model, Q, ctr, KDmodel, R_prop, T_prop, C_)
    %ICP_1 does one iteration of ICP
    
    %E step
    P_proposed = bsxfun(@plus, P*R_prop', T_prop');
    ind = knnsearch(KDmodel, P_proposed);
    P_match = P_model(ind, :);
    
    %M step
    C_model = mean(P_match, 1);
    Q_model = bsxfun(@minus, P_match, C_model);
    H = Q*Q_model;
    [U, ~, V] = svd(H);
    X = V*U';
    detX = det(X);
    
    R = X;
    T = C_model' - R*ctr;
    error = norm(P_match - P_proposed, 'fro');
    showfine
end