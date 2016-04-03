function [ hoppings ] = graphene_interlayer_hopping_Pablo1( bond_vec_21,angle_from,angle_to )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%graphene_env;

Vpppi0=-2.7;
Vppsigma0=0.48;

aa0=1.42;

ddelta=0.184*aa0;

% bond_angles=angle(bond_vec_21(:,1)+i*bond_vec_21(:,2));
% 
% theta21=bond_angles-angle_from;
% theta12=pi+bond_angles-angle_to;



u_bond(:,1:2)=bond_vec_21(:,1:2);
ss1=size(bond_vec_21);
u_bond(1:ss1(1),3)=3.35;

d_bond=sqrt(u_bond(:,1).^2+u_bond(:,2).^2+u_bond(:,3).^2);

u_bond(:,1)=u_bond(:,1)./d_bond;
u_bond(:,2)=u_bond(:,2)./d_bond;
u_bond(:,3)=u_bond(:,3)./d_bond;




hoppings=Vpppi0*exp(-(d_bond-aa0)/ddelta).*(1-u_bond(:,3).^2)+Vppsigma0*exp(-(d_bond-3.35)/ddelta).*(u_bond(:,3).^2);



end

