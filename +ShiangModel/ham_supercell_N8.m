function [ Hmat] = ham_supercell_N8( know, total_num, inplaneGAB_table_bottom,inplaneGAB_vecs_bottom, inplaneGAB_table_top,inplaneGAB_vecs_top,inplaneGAABB_table_bottom,inplaneGAABB_vecs_bottom, inplaneGAABB_table_top,inplaneGAABB_vecs_top,inplane_hops,inter_table,inter_vec,inter_hop)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

% know is 2d vector
 
% minimal model with nearest hopping.
inplane_hop1=inplane_hops(1);
inplane_hop2=inplane_hops(2);
inplane_hop3=inplane_hops(3);
inplane_hop4=inplane_hops(4);
inplane_hop5=inplane_hops(5);
inplane_hop6=inplane_hops(6);
inplane_hop7=inplane_hops(7);
inplane_hop8=inplane_hops(8);

% inplane A to B only!
% interaction, bottom to top layer only



Hmat=zeros(total_num);
Hmat_inplane=zeros(total_num);
Hmat_interaction=zeros(total_num);

% index_bottom=1:(total_num/2);
% index_top=(total_num/2+1):total_num;

% inplane part

% 3 nearest bonds and for 2 layers and 3 3rd nearest bonds

NGAB_hoppings=[inplane_hop1*[1,1,1],inplane_hop3*[1,1,1],inplane_hop4*[1,1,1,1,1,1],inplane_hop7*[1,1,1,1,1,1],inplane_hop8*[1,1,1]];


for ind=1:21
    delvec_b=inplaneGAB_vecs_bottom(ind,:)';
    delvec_t=inplaneGAB_vecs_top(ind,:)';
    
    tmp_table_b=inplaneGAB_table_bottom(:,:,ind);
    tmp_table_t=inplaneGAB_table_top(:,:,ind);
    
    
    % convert index
    tmp_index_b=total_num*(tmp_table_b(:,2)-1)+tmp_table_b(:,1);
    tmp_index_t=total_num*(tmp_table_t(:,2)-1)+tmp_table_t(:,1);

    
    Hmat_inplane(tmp_index_b)=Hmat_inplane(tmp_index_b)+NGAB_hoppings(ind)*exp(-i*dot(know,delvec_b));
    Hmat_inplane(tmp_index_t)=Hmat_inplane(tmp_index_t)+NGAB_hoppings(ind)*exp(-i*dot(know,delvec_t));
    
end
Hmat_inplane=(Hmat_inplane+Hmat_inplane');


NGAABB_hoppings=[inplane_hop2*[1,1,1,1,1,1],inplane_hop5*[1,1,1,1,1,1],inplane_hop6*[1,1,1,1,1,1]];



% inplane 2nd nearest bonds to be added!
for ind=1:18
    delvec_b=inplaneGAABB_vecs_bottom(ind,:)';
    delvec_t=inplaneGAABB_vecs_top(ind,:)';
    
    tmp_table_b=inplaneGAABB_table_bottom(:,:,ind);
    tmp_table_t=inplaneGAABB_table_top(:,:,ind);
    
    
    % convert index
    tmp_index_b=total_num*(tmp_table_b(:,2)-1)+tmp_table_b(:,1);
    tmp_index_t=total_num*(tmp_table_t(:,2)-1)+tmp_table_t(:,1);

    Hmat_inplane(tmp_index_b)=Hmat_inplane(tmp_index_b)+NGAABB_hoppings(ind)*exp(-i*dot(know,delvec_b));
    Hmat_inplane(tmp_index_t)=Hmat_inplane(tmp_index_t)+NGAABB_hoppings(ind)*exp(-i*dot(know,delvec_t));
    
end

% diagonal part? to be added later!



% interaction part

tmp_index=total_num*(inter_table(:,2)-1)+inter_table(:,1);

% Hmat_interaction(tmp_index)=Hmat_interaction(tmp_index)+inter_hop.*exp(-i*inter_vec*know');
intsize=size(inter_table);
for inds=1:intsize(1)
	Hmat_interaction(tmp_index(inds))=Hmat_interaction(tmp_index(inds))+inter_hop(inds)*exp(-i*inter_vec(inds,:)*know');
end


Hmat_interaction=Hmat_interaction+Hmat_interaction';

 Hmat=Hmat_inplane+Hmat_interaction;

% Hmat=Hmat_interaction;

end

