clear all; clc

%%%%%%%%%%%%%%%%%%%%
% preparation step %
%      M > N       %
%%%%%%%%%%%%%%%%%%%%


sc_m=8;
sc_n=6;

graphene_setup_supercell_env_Pablo1;
graphene_supercell_env_Pablo1;
build_supercell_hamN8_Pablo1; 

fprintf('Twist angle: %f\n', rot_theta/pi*180);

%%%%%%%%%%%%%%%%%%%%%%%%%%
% positions of the atoms %
%%%%%%%%%%%%%%%%%%%%%%%%%%

X=sc_all_points;
X(:,3)=0;
X=X*lattice_a;
X((tot_num/2+1):tot_num,3)=3.35;

a1=sc_t1'*lattice_a;
a2=sc_t2'*lattice_a;

% figure(1); clf
% scatter3(X(:,1),X(:,2),X(:,3));
% 
% set(gca,'FontSize',25);
% title('Where the atoms are in the supercell')
% xlabel('X (A)');
% ylabel('Y (A)');
% zlabel('Z (A)');
% axis([-inf,inf,-inf,inf,-inf,inf]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pick magnetic field and disorder strength %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

F = 0.;
W = 0.;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate the hamiltonian at gamma point without disorder %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gamma_point=[0,0];
H = ham_supercell_Pablo1( gamma_point, total_num, inplaneGAB_table_bottom,inplaneGAB_vecs_bottom, inplaneGAB_table_top,inplaneGAB_vecs_top,inplaneGAABB_table_bottom,inplaneGAABB_vecs_bottom,inplaneGAABB_table_top,inplaneGAABB_vecs_top,inplane_hops,inter_table,inter_vec,inter_hop);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Introduce onsite disorder and magnetic Peierls phase into Hamiltonian matrix %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = size(X,1);
for i=1:N
    H(i,i) = H(i,i) - 0.78 + (rand-.5)*W;
    H(i,:) = H(i,:) .* exp(1i*pi*F*(X(:,1)*X(i,2)-X(:,2)*X(i,1))).';
end
H0 = H; X0 = X;


% H = H0(1:end/2,1:end/2);
% X = X0(1:end/2,:);
figure(3); spy(H)
[V,D] = eig(H);
D = diag(D);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Phenomenological Parameters %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tau = 400.;
T = .0025;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the time-averaged current operator %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R1 = Derivation(a1, a2, X, H, 1);
R1 = V'*R1*V;

R2 = Derivation(a1, a2, X, H, 2);
R2 = V'*R2*V;

N = size(X,1);
LR1 = zeros(N);
LR2 = zeros(N);
for a = 1:N
    LR1(:,a) = R1(:,a)./(1/tau + 1i*(D(a) - D));
    LR2(:,a) = R2(:,a)./(1/tau + 1i*(D(a) - D));
end
LR1 = V*LR1*V';
LR2 = V*LR2*V';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the Fermi-Dirac ground state and the space-averaged current %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d11 = zeros(N,41);
d12 = zeros(N,41);


for i=0:80
    mu = i/10-4
    FD = V*bsxfun(@times, 1./(1+exp((D-mu)/T)), V');

    FD = Derivation(a1, a2, X, FD, 1);

    d11(:,i+1) = -sum(FD.*LR1.', 2);
    d12(:,i+1) = sum(FD.*LR2.', 2);
end

s11 = zeros(1, 81);
s12 = zeros(1, 81);
r11 = zeros(1, 81);
r12 = zeros(1, 81);
for i=1:81
   sigma = zeros(1,2);
   s11(i) = mean(d11(:,i),1);
   s12(i) = mean(d12(:,i),1);
   r11(i) = s11(i)/(s11(i)'*s11(i) + s12(i)'*s12(i))/(2*pi);
   r12(i) = s12(i)/(s11(i)'*s11(i) + s12(i)'*s12(i))/(2*pi);
end

dos = zeros(1,81);
delta = .1;
for i=0:80
    mu = i/10-4
    dos(i+1) = 1/pi*mean(imag(1./(D - mu - 1i*delta)));
end

figure(1); clf
subplot(1,3,1); hold on
plot((0:80)/10-4, real(s11))
subplot(1,3,2); hold on
plot((0:80)/10-4, real(s12))
subplot(1,3,3); hold on
plot((0:80)/10-4, sum(dos,1))

figure(4); clf
err11 = zeros(size(d11,1),1);
err12 = zeros(size(d12,1),1);

I = 42;
err11(:) = d11(:,I) - mean(d11(:,I));
err12(:) = d12(:,I) - mean(d12(:,I));
quiver3(X(:,1), X(:,2), X(:,3), real(d11(:,I)), real(d12(:,I)), zeros(size(X,1),1))
% subplot(1,2,1)
% scatter3(X(:,1), X(:,2), X(:,3), [], real(err11));
% colorbar;
% subplot(1,2,2)
% scatter3(X(:,1), X(:,2), X(:,3), [], real(err12));
% colorbar;
%histogram2(real(d12(:,20)), imag(d12(:,20)), 'FaceColor', 'flat', 'BinMethod', 'scott', 'DisplayStyle', 'tile');

