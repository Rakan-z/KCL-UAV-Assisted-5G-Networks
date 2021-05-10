function [dnet_lin, ec1_lin, ec2_lin, coordinates_optimal_edge_lin, x_optimal_edge_lin] = FBS_Boundary_search_lin(cell_diameter, GNs, UAV_parameters)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function will call an algorithm that finds the optimal linear cell edge location between two terrestrial Macro Base Stations (MBSs) based on
%the trajectory optimization of the Flying Base Stations (FBSs) which will be determined using the tspsearch function calling the 2-opt algorithm 
%- Author of tspsearch: Jonas Lundgren <splinefit@gmail.com> 2012
%----------------------------------------------------------------------------------------------------------------------------------------------------
%Inputs:
%cell_diameter - distance(km) between the two terrestrial MBSs
%GNs - #nodesX2 matrix of coordinates of ground nodes where every coordinate 0≤n≤cell_diameter
%UAV_parameters - 1X2 matrix including the UAV parameters: [power(W) velocity(m/s)]

%Outputs:
%dnet_lin - matrix containing net distance(km) covered by FBSs at every linear cell edge (boundary) location
%ec1_lin - matrix containing energy (J = Ws) consumed by FBS 1 at every linear cell edge location
%ec2_lin - matrix containing energy (J = Ws) consumed by FBS 2 at every linear cell edge location
%coordinates_optimal_edge_lin - 2X2 matrix containg coordinates of an optimal linear cell edge location
%x_optimal_edge_lin - x coordinate of an optimal linear cell edge location, used by the UAV_CellSplit algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%terrestrial MBSs 1&2 coordinates & concatenation with the matrix containing GNs to form M with all point coordinates
MBS = [0 (cell_diameter/2) ; cell_diameter (cell_diameter/2)]; 
M = [MBS;GNs];

%creating power and velocity variables using UAV_parameters input  
UAV_power = UAV_parameters(1);
UAV_velocity = UAV_parameters(2);

%creating matrix of x coordinates of all nodes
x_M = M(:,1);

%creating matrix containing inital boundary x coordinates that will be incremented with every iteration of loop  
x_boundary = [((cell_diameter/2) - (4*(cell_diameter/20))) ((cell_diameter/2) - (4*(cell_diameter/20)))];

%creating matrix containing boundary y coordinates that will remain constant
y_boundary = [0 cell_diameter];

%creating matrix containing x coordinates of all boundaries 
x_boundaries = [(x_boundary(1) + cell_diameter/20) (x_boundary(1) + (2 * cell_diameter/20)) (x_boundary(1) + (3 * cell_diameter/20)) (x_boundary(1) + (4 * cell_diameter/20)) (x_boundary(1) + (5 * cell_diameter/20)) (x_boundary(1) + (6 * cell_diameter/20)) (x_boundary(1) + (7 * cell_diameter/20))];

z=0; %initialization of figure incrementation

%loop iterating to increment boundary x coordinates to experiment with the 7 boundary locations
for B = 1:7
            
            x_boundary(1) = x_boundary(1) + cell_diameter/20; %incrementation of boundary x coordinates to iterate through different boundary locations
                        
            x_boundary(:) = x_boundary(1);
            
            X = M(x_M < x_boundary(1), :); %input matrix of tspsearch function including MBS 1 and all GNs left of the boundary (points with x < x_boundary)
            s = size(X,1); %number of nodes in X
            [p,d_left] = tspsearch(X,s); %2-opt NN algorithm with input: X,s and output: trajectory p, distance d_left of FBS 1 
            
            %sorting X by order of nodes starting from MBS 1 in trajectory p 
            p_left = X(p,:);
            p_left = [p_left;p_left(1,:)];
            
            %finding energy consumed by FBS 1
            t1 = 7+(9-7).*rand(s,1);
            hovering_t1 = sum(t1);
            ec1_lin(1,B) = UAV_power * ((d_left/UAV_velocity) + hovering_t1);

            X = M(x_M > x_boundary(1), :); %input matrix of tspsearch function including MBS 2 and all GNs right of the boundary (points with x > x_boundary)
            s = size(X,1); %number of nodes in X
            [p,d_right] = tspsearch(X,s); %2-opt NN algorithm with input: X,s and output: trajectory p, distance d_right of FBS 2
           
            %sorting X by order of nodes starting from MBS 2 in trajectory p 
            p_right = X(p,:);
            p_right = [p_right;p_right(1,:)];

            %finding energy consumed by FBS 2
            t2 = 7+(9-7).*rand(s,1);
            hovering_t2 = sum(t2);
            ec2_lin(1,B) = UAV_power * ((d_right/UAV_velocity) + hovering_t2);
            
            %creating matrix containing net distance covered by FBSs at each boundary position
            dnet_lin(1,B) = d_left + d_right;

            %plotting figure with MBSs & GNs which can be used to visualize the iterations through different boundary locations and FBS trajectories 
            z = z+1;
            figure(1) %change 1 to z to see plots at all instances
            plot(M(1,1),M(1,2),'bd',M(2,1),M(2,2),'rd',M(x_M<x_boundary(1),1),M(x_M<x_boundary(1),2),'b*',M(x_M>x_boundary(1),1),M(x_M>x_boundary(1),2),'r*',p_left(:,1),p_left(:,2),'b-',p_right(:,1),p_right(:,2),'r-',x_boundary,y_boundary,'k-','LineWidth',1)
end

%creating variables & matrices that will contain the algorithm outputs that display an optimal boundary location 
optimal_edgeIndex = find(dnet_lin==min(dnet_lin));
x_optimal_edge_lin = x_boundaries(optimal_edgeIndex(1));
coordinates_optimal_edge_lin = [x_optimal_edge_lin(1) x_optimal_edge_lin(1) ; y_boundary]';
end