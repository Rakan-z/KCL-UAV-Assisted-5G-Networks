%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%this program will run a monte carlo simulation using the FBS Boundary-search algorithm to find the most optimal linear cell edge between the two 
%Macro Base Stations based on the net distances of the two FBSs when there is an overcrowding in one cell
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nodes_i = 7; %number of iterations of node sizes (20 nodes, 30 nodes, 40 nodes ...)
montecarlo_i = 1000; %number of monte carlo instances
nodes = 20; %starting number of nodes
nodes = nodes-10; %initializing node size for the for loop incrementations 
seed1 = 0; %initializing seed 1 which will be used to maintain the distribution of points throughout the instance i of the simulation
seed2 = nodes_i*montecarlo_i + 1; %initializing seed 2 at the maximum number of iterations so that the x and y matrices below do not interchange the seeds
net_distance_lin = zeros(nodes_i*montecarlo_i, 7); %initializing matrix which will contain net distance
cell_diameter = 3; %input of distance(km) between the two MBSs
UAV_parameters = [50 60]; %input of [power(W) velocity(m/s)] of UAV

%loop 1 will increment the number of nodes to run montecarlo_i experiments for each number of nodes from 20 to 80 nodes
for node = 1:nodes_i                                                            

    nodes = nodes + 10; %incrementation of nodes

    for i = 1:montecarlo_i
        
        %incrementation of seeds so that after every instance i of the monte carlo simulation a new random set of GNs is distributed
        seed1 = seed1 + 1;
        seed2 = seed2 + 1;
        
        %fixing a probability distribution that will be used to generate 45% of GNs in a strip of area 3kmx200m representing a flashcrowd 
        cdf = 0.325; 
        pd=makedist('PiecewiseLinear','x',[0 0.65 0.85 cell_diameter],'Fx',[0 cdf*0.65 (cdf*0.85)+0.45 1]);
        
        rng(seed1) %seeding the rand generator with seed 1
        x=random(pd,nodes,1); %randomizing x coordinates of nodes GNs using the probability distribution       
        rng(seed2) %seeding the rand generator with seed 2
        y=0+(cell_diameter-0).*rand(nodes,1); %randomizing y coordinates of nodes GNs ranging from 0 to the cell diameter
        GNs = [x y];
                
        %calling function to run UAV-CellSplit_lin algorithm to find the optimal linear edge at every instance i of the monte carlo simulation
        [dnet_lin, ec1_lin, ec2_lin, coordinates_optimal_edge_lin, x_optimal_edge_lin] = FBS_Boundary_search_lin(cell_diameter, GNs, UAV_parameters);

        montecarlo_i_nodes_i = i + (montecarlo_i*(node-1)); %variable which will store the total number of runs of the simulation with every iteration  
        net_distance_lin(montecarlo_i_nodes_i,:) = dnet_lin; %creating matrix containing net distance covered by FBSs at each boundary position for every i for all nodes        
    end
end
%finding the standard deviation of net distances at each linear boundary of all instances of the monte carlo simulation for every set of nodes
dnet_std_20nodes = std(net_distance_lin(1:montecarlo_i,:));
dnet_std_30nodes = std(net_distance_lin(montecarlo_i+1:2*montecarlo_i,:));
dnet_std_40nodes = std(net_distance_lin(2*montecarlo_i+1:3*montecarlo_i,:));
dnet_std_50nodes = std(net_distance_lin(3*montecarlo_i+1:4*montecarlo_i,:));
dnet_std_60nodes = std(net_distance_lin(4*montecarlo_i+1:5*montecarlo_i,:));
dnet_std_70nodes = std(net_distance_lin(5*montecarlo_i+1:6*montecarlo_i,:));
dnet_std_80nodes = std(net_distance_lin(6*montecarlo_i+1:7*montecarlo_i,:));


%finding the mean net distance at each boundary of all instances of the monte carlo simulation for every set of nodes and the sd of the results
mean_20nodes = mean(net_distance_lin(1:montecarlo_i,:));
results_sd_20 = std(mean_20nodes);
mean_30nodes = mean(net_distance_lin(montecarlo_i+1:2*montecarlo_i,:));
results_sd_30 = std(mean_30nodes);
mean_40nodes = mean(net_distance_lin(2*montecarlo_i+1:3*montecarlo_i,:));
results_sd_40 = std(mean_40nodes);
mean_50nodes = mean(net_distance_lin(3*montecarlo_i+1:4*montecarlo_i,:));
results_sd_50 = std(mean_50nodes);
mean_60nodes = mean(net_distance_lin(4*montecarlo_i+1:5*montecarlo_i,:));
results_sd_60 = std(mean_60nodes);
mean_70nodes = mean(net_distance_lin(5*montecarlo_i+1:6*montecarlo_i,:));
results_sd_70 = std(mean_70nodes);
mean_80nodes = mean(net_distance_lin(6*montecarlo_i+1:7*montecarlo_i,:));
results_sd_80 = std(mean_80nodes);

%bar chart results for every set of nodes containing mean net distance at each boundary and highlighting the optimal boundary
figure(2)
hold on
for bc1 = 1:length(mean_20nodes)
    result_20nodes = bar(bc1,mean_20nodes(bc1),'BarWidth',0.5);
    if mean_20nodes(bc1) == min(mean_20nodes)
        set(result_20nodes,'FaceColor','g');
    else
        set(result_20nodes,'FaceColor','r');
    end
    title('20 Nodes: Optimal Linear Boundary')
    xlabel('Boundary x-coordinate (km)')
    ylabel('Mean Net-Distance (km)')
    set(gca,'xticklabel',{[],'1.05','1.2','1.35','1.5','1.65','1.8','1.95'})
end
hold off

figure(3)
hold on
for bc2 = 1:length(mean_30nodes)
    result_30nodes = bar(bc2,mean_30nodes(bc2),'BarWidth',0.5);
    if mean_30nodes(bc2) == min(mean_30nodes)
        set(result_30nodes,'FaceColor','g');
    else
        set(result_30nodes,'FaceColor','r');
    end
    title('30 Nodes: Optimal Linear Boundary')
    xlabel('Boundary x-coordinate (km)')
    ylabel('Mean Net-Distance (km)')
    set(gca,'xticklabel',{[],'1.05','1.2','1.35','1.5','1.65','1.8','1.95'})
end
hold off 

figure(4)
hold on
for bc3 = 1:length(mean_40nodes)
    result_40nodes = bar(bc3,mean_40nodes(bc3),'BarWidth',0.5);
    if mean_40nodes(bc3) == min(mean_40nodes)
        set(result_40nodes,'FaceColor','g');
    else
        set(result_40nodes,'FaceColor','r');
    end
    title('40 Nodes: Optimal Linear Boundary')
    xlabel('Boundary x-coordinate (km)')
    ylabel('Mean Net-Distance (km)')
    set(gca,'xticklabel',{[],'1.05','1.2','1.35','1.5','1.65','1.8','1.95'})
end
hold off 

figure(5)
hold on
for bc4 = 1:length(mean_50nodes)
    result_50nodes = bar(bc4,mean_50nodes(bc4),'BarWidth',0.5);
    if mean_50nodes(bc4) == min(mean_50nodes)
        set(result_50nodes,'FaceColor','g');
    else
        set(result_50nodes,'FaceColor','r');
    end
    title('50 Nodes: Optimal Linear Boundary')
    xlabel('Boundary x-coordinate (km)')
    ylabel('Mean Net-Distance (km)')
    set(gca,'xticklabel',{[],'1.05','1.2','1.35','1.5','1.65','1.8','1.95'})
end
hold off

figure(6)
hold on
for bc5 = 1:length(mean_60nodes)
    result_60nodes = bar(bc5,mean_60nodes(bc5),'BarWidth',0.5);
    if mean_60nodes(bc5) == min(mean_60nodes)
        set(result_60nodes,'FaceColor','g');
    else
        set(result_60nodes,'FaceColor','r');
    end
    title('60 Nodes: Optimal Linear Boundary')
    xlabel('Boundary x-coordinate (km)')
    ylabel('Mean Net-Distance (km)')
    set(gca,'xticklabel',{[],'1.05','1.2','1.35','1.5','1.65','1.8','1.95'})
end
hold off

figure(7)
hold on
for bc6 = 1:length(mean_70nodes)
    result_70nodes = bar(bc6,mean_70nodes(bc6),'BarWidth',0.5);
    if mean_70nodes(bc6) == min(mean_70nodes)
        set(result_70nodes,'FaceColor','g');
    else
        set(result_70nodes,'FaceColor','r');
    end
    title('70 Nodes: Optimal Linear Boundary')
    xlabel('Boundary x-coordinate (km)')
    ylabel('Mean Net-Distance (km)')
    set(gca,'xticklabel',{[],'1.05','1.2','1.35','1.5','1.65','1.8','1.95'})
end
hold off

figure(8)
hold on
for bc7 = 1:length(mean_80nodes)
    result_80nodes = bar(bc7,mean_80nodes(bc7),'BarWidth',0.5);
    if mean_80nodes(bc7) == min(mean_80nodes)
        set(result_80nodes,'FaceColor','g');
    else
        set(result_80nodes,'FaceColor','r');
    end
    title('80 Nodes: Optimal Linear Boundary')
    xlabel('Boundary x-coordinate (km)')
    ylabel('Mean Net-Distance (km)')
    set(gca,'xticklabel',{[],'1.05','1.2','1.35','1.5','1.65','1.8','1.95'})
end
hold off