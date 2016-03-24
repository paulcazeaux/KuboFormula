% Number of points %
N = 80;

% Build grid points %
[X,Y] = meshgrid((1:N) - floor(.5*(N+1)));
X = [X(:), Y(:)];
clear Y;

% Connectivity %
C = zeros(N*N, 4);
for i=1:N*N
    [~, C(i,1)] = ismember(mod(X(i,:)+[1,0],N), mod(X,N), 'rows');
    [~, C(i,2)] = ismember(mod(X(i,:)+[0,1],N), mod(X,N), 'rows');
    [~, C(i,3)] = ismember(mod(X(i,:)+[-1,0],N), mod(X,N), 'rows');
    [~, C(i,4)] = ismember(mod(X(i,:)+[0,-1],N), mod(X,N), 'rows');
end

%% Pick magnetic field and disorder strength %%
F = 0*2/N;
W = 0;

% Build Hamiltonian matrix %
H = zeros(N*N);
for i=1:N*N
    H(i,i) = (rand-.5)*W;
    H(i,C(i,1)) = 1*exp(1i*pi*F*(1 *X(i,2) -0*X(i,1)));
    H(i,C(i,2)) = 1*exp(1i*pi*F*(0 *X(i,2) -1*X(i,1)));
    H(i,C(i,3)) = 1*exp(1i*pi*F*(-1*X(i,2) -0*X(i,1)));
    H(i,C(i,4)) = 1*exp(1i*pi*F*(0 *X(i,2) +1*X(i,1)));
end

[V,D] = eig(H);
D = diag(D);

R1 = Derivation(N, X, H, 1);
R1 = V'*R1*V;

R2 = Derivation(N, X, H, 2);
R2 = V'*R2*V;

%% Phenomenological Parameters %%
mu = -4/9;
T = .1;
tau = 10.;

FD = V*diag(1./(1+exp((D-mu)/T)))*V';
S1 = Derivation(N, X, FD, 1);
S1 = V'*S1*V;

S2 = Derivation(N, X, FD, 2);
S2 = V'*S2*V;

sigma = zeros(2);

Ntot = N*N;
for a = 1:Ntot
    for b = 1:Ntot
        sigma(1,1) = sigma(1,1) + R1(b,a)*S1(a,b)/(1/tau + 1i*(D(a) - D(b)));
        sigma(1,2) = sigma(1,2) + R1(b,a)*S2(a,b)/(1/tau + 1i*(D(a) - D(b)));
        sigma(2,1) = sigma(2,1) + R2(b,a)*S1(a,b)/(1/tau + 1i*(D(a) - D(b)));
        sigma(2,2) = sigma(2,2) + R2(b,a)*S2(a,b)/(1/tau + 1i*(D(a) - D(b)));
    end
end
sigma = -sigma/Ntot;
disp(sigma)