function res = Derivation(a1, a2, X, G, j)
    % First write the beta coefficients, solution of 
    % a Vandermonde linear system
    Q = 10; % Assume that N > 2Q+1
    beta = zeros(Q,1);
    for i=1:Q
        beta(i) = i*prod(1 - (i./[1:i-1 i+1:Q]).^2);
    end
    
    
    beta = 1./(4*pi*beta);
    N = X(:, 1:2)/[a1; a2];
    P1 = exp(2*1i*pi*N(:,1)); % P1, P2 are column vectors
    P2 = exp(2*1i*pi*N(:,2));
    
    res = zeros(size(G));
    for k=1:length(beta)
        M1 = bsxfun(@times, bsxfun(@times, P1.^(-k),  G),   P1.^k.' );
        M2 = bsxfun(@times, bsxfun(@times, P1.^k,     G), P1.^(-k).');
        res = res + a1(j)*beta(k)*(M1 - M2);
        
        M1 = bsxfun(@times, bsxfun(@times, P2.^(-k),  G),   P2.^k.' );
        M2 = bsxfun(@times, bsxfun(@times, P2.^k,     G), P2.^(-k).');
        res = res + a2(j)*beta(k)*(M1 - M2);
    end
end