function res = Derivation(N, X, G, j)
    % First write the beta coefficients, solution of 
    % a Vandermonde linear system
    Q = 10; % Assume that N > 2Q+1
    beta = zeros(Q,1);
    for i=1:Q
        beta(i) = i*prod(1 - (i./[1:i-1 i+1:Q]).^2);
    end
    beta = N./(4*pi*beta);
    
    P = exp(2*1i*pi*X(:,j)/N); % P is a column vector
    
    res = zeros(size(G));
    for k=1:length(beta)
        M1 = bsxfun(@times, bsxfun(@times, P.^k,     G), P.^(-k).');
        M2 = bsxfun(@times, bsxfun(@times, P.^(-k),  G),   P.^k.' );
        res = res + beta(k)*(M1 - M2);
    end
end