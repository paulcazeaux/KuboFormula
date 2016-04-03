clear all;


% M > N
% preparation step
% M > N
 sc_m=4*1;
 sc_n=3*1; 
% 
graphene_setup_supercell_env
graphene_supercell_env

% figure(1);
% supercell_plot;
build_supercell_hamN8;

'twist angle'

rot_theta/pi*180

% positions of the atoms
pos_all_points=sc_all_points;


pos_all_points(:,3)=0;
pos_all_points=pos_all_points*lattice_a;

pos_all_points((tot_num/2+1):tot_num,3)=3.35;

pos_a1=sc_t1'*lattice_a;
pos_a2=sc_t2'*lattice_a;

scatter3(pos_all_points(:,1),pos_all_points(:,2),pos_all_points(:,3));
title('Where the atoms are in the supercell')
xlabel('X (A)');
ylabel('Y (A)');
zlabel('Z (A)');
set(gca,'FontSize',25);
axis([-inf,inf,-inf,inf,-inf,inf]);

% generate the hamiltonian at gamma point without disorder
gamma_point=[0,0];
Hmat_Gamma = ham_supercell_N8( gamma_point, total_num, inplaneGAB_table_bottom,inplaneGAB_vecs_bottom, inplaneGAB_table_top,inplaneGAB_vecs_top,inplaneGAABB_table_bottom,inplaneGAABB_vecs_bottom,inplaneGAABB_table_top,inplaneGAABB_vecs_top,inplane_hops,inter_table,inter_vec,inter_hop);


% where Paul can start with conductance computation!



