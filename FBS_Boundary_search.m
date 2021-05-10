function [dnet_pw, ec1_pw, ec2_pw, coordinates_optimal_edge_pw] = FBS_Boundary_search(cell_diameter, GNs, UAV_parameters)

%This function will call an algorithm that finds the optimal piecewise cell edge location between two terrestrial Macro Base Stations (MBSs) based on
%the location of the optimal linear cell edge and the trajectory optimization of the Flying Base Stations (FBSs) which will be determined using the 
%tspsearch function calling the 2-opt algorithm - Author of tspsearch: Jonas Lundgren <splinefit@gmail.com> 2012

%Inputs:
%cell_diameter - distance(km) between the two terrestrial MBSs
%GNs - #nodesX2 matrix of coordinates of ground nodes where every coordinate 0≤n≤cell_diameter
%UAV_parameters - 1X2 matrix including the UAV parameters: [power(W) velocity(km/h)]

%Outputs:
%dnet_pw - matrix containing net distance(km) covered by FBSs at every piecewise cell edge location
%ec1_pw - matrix containing energy (J = Ws) consumed by FBS 1 at every piecewise cell edge location
%ec2_pw - matrix containing energy (J = Ws) consumed by FBS 2 at every piecewise cell edge location
%coordinates_optimal_edge_pw - 2X2 matrix containg coordinates of an optimal piecewise cell edge location


[dnet_lin, ec1_lin, ec2_lin, coordinates_optimal_edge_lin, x_optimal_edge_lin] = FBS_Boundary_search_lin(cell_diameter, GNs, UAV_parameters)
[dnet_pw, ec1_pw, ec2_pw, coordinates_optimal_edge_pw] = FBS_Boundary_search_pw(cell_diameter, GNs, UAV_parameters, x_optimal_edge_lin)
end
