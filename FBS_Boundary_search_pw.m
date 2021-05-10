function [dnet_pw, ec1_pw, ec2_pw, coordinates_optimal_edge_pw] = FBS_Boundary_search_pw(cell_diameter, GNs, UAV_parameters, x_optimal_edge_lin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function will call an algorithm that finds the optimal piecewise cell edge location between two terrestrial Macro Base Stations (MBSs) based on
%the location of the optimal linear cell edge and the trajectory optimization of the Flying Base Stations (FBSs) which will be determined using the 
%tspsearch function calling the 2-opt algorithm - Author of tspsearch: Jonas Lundgren <splinefit@gmail.com> 2012
%----------------------------------------------------------------------------------------------------------------------------------------------------
%Inputs:
%cell_diameter - distance(km) between the two terrestrial MBSs
%GNs - #nodesX2 matrix of coordinates of ground nodes where every coordinate 0≤n≤cell_diameter
%UAV_parameters - 1X2 matrix including the UAV parameters: [power(W) velocity(m/s)]
%x_optimal_edge_lin - x coordinate of the optimal linear cell edge (boundary) used to find the range of points which the piecewise cell edge will 
%iterate between 

%Outputs:
%dnet_pw - matrix containing net distance(km) covered by FBSs at every piecewise cell edge location
%ec1_pw - matrix containing energy (J = Ws) consumed by FBS 1 at every piecewise cell edge location
%ec2_pw - matrix containing energy (J = Ws) consumed by FBS 2 at every piecewise cell edge location
%coordinates_optimal_edge_pw - 2X2 matrix containg coordinates of an optimal piecewise cell edge location
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%terrestrial MBSs 1&2 coordinates & concatenation with the matrix containing GNs to form M with all point coordinates
MBS = [0 (cell_diameter/2) ; cell_diameter (cell_diameter/2)]; 
M = [MBS;GNs];

%creating power and velocity variables using UAV_parameters input  
UAV_power = UAV_parameters(1);
UAV_velocity = UAV_parameters(2);

%number of locations that two of the piecewise boundary points (p2,p3) will iterate through to find trajectories at multiple piecewise boundaries
p2_i = 5;
p3_i = 5;

%to guarantee that the optimal piecewise boundaries will be found for one of the optimal linear boundaries incase there was more than one
if length(x_optimal_edge_lin) > 1
    x_optimal_edge_lin = x_optimal_edge_lin(1);
end

%splitting M into two matrices representing x and y coordinates of all points 
x = M(:,1);
y = M(:,2);

%creating matrix containing inital boundary x coordinates that will be incremented with every iteration of loop  
x_boundary_pw = [x_optimal_edge_lin (x_optimal_edge_lin - (3 * cell_diameter/40)) (x_optimal_edge_lin - (3 * cell_diameter/40)) x_optimal_edge_lin];

%creating matrix containing boundary y coordinates that will remain constant
y_boundary_pw = [cell_diameter (cell_diameter/4 + cell_diameter/2) (cell_diameter/2 - cell_diameter/4) 0];

%creating matrix containing x coordinates of all locations of points that piecewise boundary will iterate through 
xP2P3_cell_edges_pw = [(x_boundary_pw(2) + cell_diameter/40) (x_boundary_pw(2) + (2 * cell_diameter/40)) (x_boundary_pw(2) + (3 * cell_diameter/40)) (x_boundary_pw(2) + (4 * cell_diameter/40)) (x_boundary_pw(2) + (5 * cell_diameter/40))];

u=0; %initialization of figure incrementation


%loops iterating to increment boundary (p2,p3) x coordinates to experiment with the 25 boundary locations
%at every p2 there will be 5 incrementations of p3 
for p2 = 1:p2_i
            
            %incrementation of boundary p2 x coordinate to iterate through different boundary locations
            x_boundary_pw(2) = x_boundary_pw(2) + cell_diameter/40; 
            
            %initializing p3 at first location after the second loop increments p3 to the upper bound of the range of the x coordinates of the boundary
            x_boundary_pw(3) = x_boundary_pw(1) - (3*cell_diameter/40); 
            
            
            for p3 = 1:p3_i

                p2_i_p3_i = p3 + (p3_i*(p2-1)); %variable which will store the position of the boundary from 1-25 with every iteration 
                
                %incrementation of boundary p3 x coordinate to iterate through different boundary locations            
                x_boundary_pw(3) = x_boundary_pw(3) + cell_diameter/40;

                %building polygon which wraps around polygonal area to the left of the
                %boundary
                x_Poly_left_b(2:3) = x_boundary_pw(2:3);     
                x_Poly_left_b(1) = x_boundary_pw(1);
                x_Poly_left_b(4) = x_boundary_pw(4);
                x_Poly_left_b(5:6) = 0;
                y_Poly_left_b(6) = y_boundary_pw(1);
                y_Poly_left_b(5) = 0;
                y_Poly_left_b(1:4) = y_boundary_pw(:);
                
                %implementation of inpolygon function which determines what points of a matrix are inside a polygon to determine which points are to the left of the boundary
                left = inpolygon(x,y,x_Poly_left_b,y_Poly_left_b);

                %creating matrices containing coordinates of points on each side of the boundary
                x_nodes_left_b = x(left); %x coordinates of points to the left of the boundary 
                y_nodes_left_b = y(left); %y coordinates of points to the left of the boundary 
                s_nodes_left_b = size(x_nodes_left_b,1); %number of points to the left of the boundary

                x_nodes_right_b = x(~left); %x coordinates of points to the right of the boundary        
                y_nodes_right_b = y(~left); %y coordinates of points to the right of the boundary
                s_nodes_right_b = size(x_nodes_right_b,1); %number of points to the right of the boundary

                
                X = [x_nodes_left_b y_nodes_left_b]; %input matrix of tspsearch function including MBS 1 and all GNs left of the boundary (points inside polygon)
                s = s_nodes_left_b; %number of nodes in X
                [p,d_left] = tspsearch(X,s); %2-opt NN algorithm with input: X,s and output: trajectory p, distance d_left of FBS 1
                
                %sorting X by order of nodes starting from MBS 1 in trajectory p 
                p_left = X(p,:);
                p_left = [p_left;p_left(1,:)];
                
                %finding energy consumed by FBS 1 in J =(Ws)
                t1 = 7+(9-7).*rand(s,1);
                hovering_t1 = sum(t1);
                ec1_pw(1,p2_i_p3_i) = UAV_power * ((d_left/UAV_velocity) + hovering_t1);
 
                X = [x_nodes_right_b y_nodes_right_b]; %input matrix of tspsearch function including MBS 2 and all GNs right of the boundary (points outside polygon)
                s = s_nodes_right_b; %number of nodes in X
                [p,d_right] = tspsearch(X,s); %2-opt NN algorithm with input: X,s and output: trajectory p, distance d_right of FBS 2
                
                %sorting X by order of nodes starting from MBS 2 in trajectory p 
                p_right = X(p,:);
                p_right = [p_right;p_right(1,:)];
                
                %finding energy consumed by FBS 2 in J =(Ws)
                t2 = 7+(9-7).*rand(s,1);
                hovering_t2 = sum(t2);
                ec2_pw(1,p2_i_p3_i) = UAV_power * ((d_right/UAV_velocity) + hovering_t2);
            
                dnet_pw(1,p2_i_p3_i)= d_left + d_right; %creating matrix containing net distance covered by FBSs at each boundary position

                %plotting figure with MBSs & GNs which can be used to visualize the iterations through different boundary locations and FBS trajectories 
                u = u+1;
                figure(1) %change 2 to u to see plots at all instances 
                plot(x_nodes_left_b(1),y_nodes_left_b(1),'bd',x_nodes_right_b(1),y_nodes_right_b(1),'rd',x_nodes_left_b(2:s_nodes_left_b),y_nodes_left_b(2:s_nodes_left_b),'b*',x_nodes_right_b(2:s_nodes_right_b),y_nodes_right_b(2:s_nodes_right_b),'r*',p_left(:,1),p_left(:,2),'b-',p_right(:,1),p_right(:,2),'r-',x_boundary_pw,y_boundary_pw,'k-','LineWidth',1);

            end
end

%creating matrix that will contain the algorithm output that displays an optimal boundary location 
optimal_edgeIndex_pw = find(dnet_pw==min(dnet_pw));
optimal_edgeIndex_pw = optimal_edgeIndex_pw(1);

if any(mod(optimal_edgeIndex_pw,5) == 0)
    f = find(mod(optimal_edgeIndex_pw,5) == 0);
    xP2_optimal_edge_pw = xP2P3_cell_edges_pw(optimal_edgeIndex_pw(f)/5);
    xP3_optimal_edge_pw = xP2P3_cell_edges_pw(5);
else xP2_optimal_edge_pw = xP2P3_cell_edges_pw((optimal_edgeIndex_pw - mod(optimal_edgeIndex_pw,5))/5 + 1);
     xP3_optimal_edge_pw = xP2P3_cell_edges_pw(mod(optimal_edgeIndex_pw,5));
end

coordinates_optimal_edge_pw = [x_optimal_edge_lin xP2_optimal_edge_pw xP3_optimal_edge_pw x_optimal_edge_lin;y_boundary_pw]';

end
