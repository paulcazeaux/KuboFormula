% Number of points %
N = 100;

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

% Pick magnetic field and disorder strength %%
F = 0.1;
W = 1;

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

% Taking the trace in real space
tau = 40.;

LR1 = zeros(N*N);
LR2 = zeros(N*N);
for a = 1:N*N
    LR1(:,a) = R1(:,a)./(1/tau + 1i*(D(a) - D));
    LR2(:,a) = R2(:,a)./(1/tau + 1i*(D(a) - D));
end
LR1 = V*LR1*V';
LR2 = V*LR2*V';

%% Phenomenological Parameters %%

T = .025;

d11 = zeros(N*N,41);
d12 = zeros(N*N,41);
for i=0:40
    mu = i/10-4
    FD = V*bsxfun(@times, 1./(1+exp((D-mu)/T)), V');

    FD = Derivation(N, X, FD, 1);

    d11(:,i+1) = -sum(FD.*LR1.', 2);
    d12(:,i+1) = sum(FD.*LR2.', 2);
end

%%
r11 = zeros(1, 41);
r12 = zeros(1, 41);
for i=1:41
   sigma = zeros(1,2);
   sigma(1) = mean(d11(:,i),1);
   sigma(2) = mean(d12(:,i),1);
   r11(i) = sigma(1)/(sigma(1)'*sigma(1) + sigma(2)'*sigma(2))/(2*pi);
   r12(i) = sigma(2)/(sigma(1)'*sigma(1) + sigma(2)'*sigma(2))/(2*pi);
end

figure(1); clf
subplot(1,2,1)
histogram2(real(d11(:,20)), 2*pi*imag(d11(:,20)), 'FaceColor', 'flat', 'BinMethod', 'scott', 'DisplayStyle', 'tile');
subplot(1,2,2)
histogram2(real(d12(:,20)), 2*pi*imag(d12(:,20)), 'FaceColor', 'flat', 'BinMethod', 'scott', 'DisplayStyle', 'tile');
figure(2); clf
subplot(1,2,1)
plot((6:40)/10-4, real(r11(7:41)))
subplot(1,2,2)
plot((6:40)/10-4, real(r12(7:41)))
