function [ hoppings ] = graphene_interlayer_hopping( bond_vec_21,angle_from,angle_to )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
 
graphene_env;

bond_angles=angle(bond_vec_21(:,1)+i*bond_vec_21(:,2));

theta21=bond_angles-angle_from;
theta12=pi+bond_angles-angle_to;

bond_distance=sqrt(sum(bond_vec_21.^2,2));


hoppings=V0_func_pzonly(bond_distance) + V3_func_pzonly(bond_distance).*(cos(3*theta21)+cos(3*theta12))+ V6_func_pzonly(bond_distance).*(cos(6*theta21)+cos(6*theta12));

end

