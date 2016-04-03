% setup basic environment for supercell

% input: sc_m, sc_n

% expect sc_m >=2    sc_n>=1


lattice_a=1.42*sqrt(3);

layer_d=[0,0,3.35];


a1=[sqrt(3); -1/2];
a2=[sqrt(3); +1/2];

sc_int_bottom_a1=[sc_n,sc_m];
sc_int_bottom_a2=[-sc_m,sc_n+sc_m];

sc_int_top_a1=[sc_m,sc_n];
sc_int_top_a2=[-sc_n,sc_n+sc_m];


num_pc=(sc_m^2+sc_m*sc_n+sc_n^2);
total_num=4*num_pc;

rot_theta=acos((sc_m^2+sc_n^2+4*sc_m*sc_n)/2/(sc_m^2+sc_m*sc_n+sc_n^2));
rot_mat=[cos(rot_theta),-sin(rot_theta);sin(rot_theta),cos(rot_theta)];

ra1=rot_mat*a1;
ra2=rot_mat*a2;

ra_mat=zeros(2);
ra_mat(:,1)=ra1;
ra_mat(:,2)=ra2;


% for both bottom and top unit
sc_t1=sc_n*a1+sc_m*a2;
sc_t2=-sc_m*a1+(sc_m+sc_n)*a2;


sc_ft1=[sc_t1(1),sc_t1(2),0];
sc_ft2=[sc_t2(1),sc_t2(2),0];
sc_ft3=[0,0,1];

sc_v=abs(dot(sc_ft1,cross(sc_ft2,sc_ft3)));
sc_b1=2*pi*cross(sc_ft2,sc_ft3)/sc_v;
sc_b2=2*pi*cross(sc_ft3,sc_ft1)/sc_v;
sc_b3=2*pi*cross(sc_ft1,sc_ft2)/sc_v;

sc_vec1=sc_b1(1:2);
sc_vec2=sc_b2(1:2)+sc_b1(1:2);
sc_gamma=sc_vec1*0;
sc_kpoint=(sc_vec1+sc_vec2)/3;
sc_mpoint=sc_vec1*0.5;
