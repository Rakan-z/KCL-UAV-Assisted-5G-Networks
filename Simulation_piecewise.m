%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%this program will run a monte carlo simulation using the FBS Boundary-search algorithm to find the most optimal piecewise cell edge between the 
%two Macro Base Stations based on the location of the most optimal linear cell edge using the net distances of the two FBSs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


nodes_i = 7;  %number of iterations of node sizes (20 nodes, 30 nodes, 40 nodes ...)
montecarlo_i = 1000; %number of monte carlo instances
nodes = 20; %starting number of nodes
nodes = nodes-10; %initializing node size for the for loop incrementations 
seed = 0; %initializing seed which will be used to maintain the distribution of points throughout the instance i of the simulation
net_distance_pw = zeros(nodes_i*montecarlo_i, 25); %initializing matrix which will contain net distance 
cell_diameter = 3; %input of distance(km) between the two MBSs
UAV_parameters = [50 60]; %input of [power(W) velocity(m/s)] of UAV
UAV_power = UAV_parameters(1); %power
UAV_velocity = UAV_parameters(2); %velocity
x_optimal_edge_lin = 1.5; %output of FBS_Boundary_search_lin: x coordinates of optimal linear boundary
p = 0; %initialization of figure incrementation

%loop 1 will increment the number of nodes to run montecarlo_i experiments for each number of nodes from 20 to 80 nodes
for node = 1:nodes_i                                                            

    nodes = nodes + 10; %incrementation of nodes

    for i = 1:montecarlo_i

        seed = seed + 1; %incrementation of seed so that after every instance i of the monte carlo simulation a new random set of GNs is distributed

        rng(seed); %seeding the rand generator 
        GNs = 0+(cell_diameter-0).*rand(nodes,2); %randomizing nodes GNs ranging from 0 to the cell diameter
        
        %calling function to run FBS_Boundary_search_pw algorithm to find the optimal piecewise edge at every instance i of the monte carlo simulation
        [dnet_pw, ec1_pw, ec2_pw, coordinates_optimal_edge_pw] = FBS_Boundary_search_pw(cell_diameter, GNs, UAV_parameters, x_optimal_edge_lin);


        montecarlo_i_nodes_i = i + (montecarlo_i*(node-1)); %variable which will store the total number of runs of the simulation with every iteration 
        net_distance_pw(montecarlo_i_nodes_i,:) = dnet_pw;  %creating matrix containing net distance covered by FBSs at each boundary position for every i for all nodes 
        energy_consumption1_pw(montecarlo_i_nodes_i,:) = ec1_pw; %creating matrix containing energy consumed by FBS 1 at each boundary position for every i for all nodes
        energy_consumption2_pw(montecarlo_i_nodes_i,:) = ec2_pw;%creating matrix containing energy consumed by FBS 2 at each boundary position for every i for all nodes
        
        %kmeans algorithm which will replicate 5 times to split Kmeans_GNs into 2 clusters representing the areas covered by each FBS
        opts = statset('Display','final');
        [cidx, ctrs] = kmeans(GNs, 2, 'Distance','city','Replicates',5, 'Options',opts);
                              
        %terrestrial MBSs 1&2 coordinates 
        MBS = [0 cell_diameter/2;cell_diameter cell_diameter/2];
        
        %if statement which will swap MBS matrix rows so that FBS trajectories do not cross in the case that MBS 1 is closer to centroid of cluster 2 than to centroid of cluster 1
        MBS_ctrs1 = [MBS(1,:);ctrs(1,:)];
        MBS_ctrs2 = [MBS(1,:);ctrs(2,:)];
        if pdist(MBS_ctrs1,'euclidean') > pdist(MBS_ctrs2,'euclidean')
            MBS = [cell_diameter cell_diameter/2;0 cell_diameter/2];
        else
        end
        
        X = [MBS(1,:);GNs(cidx==1,:)]; %input matrix of tspsearch function including MBS 1 and all GNs in cluster 1
        s = size(X,1); %number of nodes in X
        [p,d_left] = tspsearch(X,s); %2-opt algorithm with input: X,s and output: trajectory p, distance d_left of FBS 1
        
        %sorting X by order of nodes starting from MBS 1 in trajectory p 
        p_left = X(p,:);;
        p_left = [p_left;p_left(1,:)];
        
        %finding energy consumed by FBS 1 in J =(Ws)
        t1 = 7+(9-7).*rand(s,1);
        hovering_t1_kmeans = sum(t1);
        ec1_kmeans(montecarlo_i_nodes_i,1) = UAV_power * ((d_left/UAV_velocity) + hovering_t1_kmeans); %creating matrix containing energy consumed by FBS 1 at each clustering for every i for all nodes

        X = [MBS(2,:);GNs(cidx==2,:)]; %input matrix of tspsearch function including MBS 2 and all GNs in cluster 2
        s = size(X,1); %number of nodes in X
        [p,d_right] = tspsearch(X,s); %2-opt algorithm with input: X,s and output: trajectory p, distance d_right of FBS 2
       
        %sorting X by order of nodes starting from MBS 1 in trajectory p 
        p_right = X(p,:);
        p_right = [p_right;p_right(1,:)];
        
        %finding energy consumed by FBS 2 in J =(Ws)
        t2 = 7+(9-7).*rand(s,1);
        hovering_t2_kmeans = sum(t2);
        ec2_kmeans(montecarlo_i_nodes_i,1) = UAV_power * ((d_right/UAV_velocity) + hovering_t2_kmeans); %creating matrix containing energy consumed by FBS 2 at each clustering for every i for all nodes

        dKmeans(montecarlo_i_nodes_i,1) = d_left + d_right;  %creating matrix containing net distance covered by FBSs at each boundary position between the two clusters

        %plotting figure with MBSs & GNs which can be used to visualize the iterations through different boundary locations and FBS trajectories displaying the K-means clustering output
        p = p+1 
        figure(3) %change 3 to p to see plots at all instances
        plot(MBS(1,1),MBS(1,2),'bd',MBS(2,1),MBS(2,2),'rd',GNs(cidx==1,1),GNs(cidx==1,2),'b*',GNs(cidx==2,1),GNs(cidx==2,2),'r*',p_left(:,1),p_left(:,2),'b-',p_right(:,1),p_right(:,2),'r-','LineWidth',1)
    end
end

%visualization of distance data---------------------------------------------------------------------------------------------------------------------------------------

%Finding the standard deviation of net distances at each piecewise boundary of all instances of the monte carlo simulation for every set of nodes
dnet_std_20nodes = std(net_distance_pw(1:montecarlo_i,:));
dnet_std_30nodes = std(net_distance_pw(montecarlo_i+1:2*montecarlo_i,:));
dnet_std_40nodes = std(net_distance_pw(2*montecarlo_i+1:3*montecarlo_i,:));
dnet_std_50nodes = std(net_distance_pw(3*montecarlo_i+1:4*montecarlo_i,:));
dnet_std_60nodes = std(net_distance_pw(4*montecarlo_i+1:5*montecarlo_i,:));
dnet_std_70nodes = std(net_distance_pw(5*montecarlo_i+1:6*montecarlo_i,:));
dnet_std_80nodes = std(net_distance_pw(6*montecarlo_i+1:7*montecarlo_i,:));


%finding the standard deviation of net distances at each boundary caused by K-means clustering of all instances of the monte carlo simulation for every set of nodes
dnet_std_20nodes_kmeans = std(dKmeans(1:montecarlo_i,:));
dnet_std_30nodes_kmeans = std(dKmeans(montecarlo_i+1:2*montecarlo_i,:));
dnet_std_40nodes_kmeans = std(dKmeans(2*montecarlo_i+1:3*montecarlo_i,:));
dnet_std_50nodes_kmeans = std(dKmeans(3*montecarlo_i+1:4*montecarlo_i,:));
dnet_std_60nodes_kmeans = std(dKmeans(4*montecarlo_i+1:5*montecarlo_i,:));
dnet_std_70nodes_kmeans = std(dKmeans(5*montecarlo_i+1:6*montecarlo_i,:));
dnet_std_80nodes_kmeans = std(dKmeans(6*montecarlo_i+1:7*montecarlo_i,:));


%finding the mean net distance at each piecewise boundary of all instances of the monte carlo simulation for every set of nodes and the sd of the results
dnet_mean_20nodes = mean(net_distance_pw(1:montecarlo_i,:));
results_sd_20 = std(dnet_mean_20nodes);
dnet_mean_30nodes = mean(net_distance_pw(montecarlo_i+1:2*montecarlo_i,:));
results_sd_30 = std(dnet_mean_30nodes);
dnet_mean_40nodes = mean(net_distance_pw(2*montecarlo_i+1:3*montecarlo_i,:));
results_sd_40 = std(dnet_mean_40nodes);
dnet_mean_50nodes = mean(net_distance_pw(3*montecarlo_i+1:4*montecarlo_i,:));
results_sd_50 = std(dnet_mean_50nodes);
dnet_mean_60nodes = mean(net_distance_pw(4*montecarlo_i+1:5*montecarlo_i,:));
results_sd_60 = std(dnet_mean_60nodes);
dnet_mean_70nodes = mean(net_distance_pw(5*montecarlo_i+1:6*montecarlo_i,:));
results_sd_70 = std(dnet_mean_70nodes);
dnet_mean_80nodes = mean(net_distance_pw(6*montecarlo_i+1:7*montecarlo_i,:));
results_sd_80 = std(dnet_mean_80nodes);

%finding the mean net distance at each boundary caused by K-means clustering of all instances of the monte carlo simulation for every set of nodes
dnet_mean_20nodes_Kmeans = mean(dKmeans(1:montecarlo_i,:));
dnet_mean_30nodes_Kmeans = mean(dKmeans(montecarlo_i+1:2*montecarlo_i,:));
dnet_mean_40nodes_Kmeans = mean(dKmeans(2*montecarlo_i+1:3*montecarlo_i,:));
dnet_mean_50nodes_Kmeans = mean(dKmeans(3*montecarlo_i+1:4*montecarlo_i,:));
dnet_mean_60nodes_Kmeans = mean(dKmeans(4*montecarlo_i+1:5*montecarlo_i,:));
dnet_mean_70nodes_Kmeans = mean(dKmeans(5*montecarlo_i+1:6*montecarlo_i,:));
dnet_mean_80nodes_Kmeans = mean(dKmeans(6*montecarlo_i+1:7*montecarlo_i,:));

%visualization of energy consumed by FBS 1 data---------------------------------------------------------------------------------------------------------------------------------------

%finding the mean energy consumed by FBS 1 at each piecewise boundary of all instances of the monte carlo simulation for every set of nodes
ec1_mean_20nodes = mean(energy_consumption1_pw(1:montecarlo_i,:));
ec1_mean_30nodes = mean(energy_consumption1_pw(montecarlo_i+1:2*montecarlo_i,:));
ec1_mean_40nodes = mean(energy_consumption1_pw(2*montecarlo_i+1:3*montecarlo_i,:));
ec1_mean_50nodes = mean(energy_consumption1_pw(3*montecarlo_i+1:4*montecarlo_i,:));
ec1_mean_60nodes = mean(energy_consumption1_pw(4*montecarlo_i+1:5*montecarlo_i,:));
ec1_mean_70nodes = mean(energy_consumption1_pw(5*montecarlo_i+1:6*montecarlo_i,:));
ec1_mean_80nodes = mean(energy_consumption1_pw(6*montecarlo_i+1:7*montecarlo_i,:));

%creating matrix containing mean energy consumed by FBS 1 at the optimal boundary for each set of nodes 
ec1_opt_allnodes = [min(ec1_mean_20nodes) min(ec1_mean_30nodes) min(ec1_mean_40nodes) min(ec1_mean_50nodes) min(ec1_mean_60nodes) min(ec1_mean_70nodes) min(ec1_mean_80nodes)];

%finding the mean energy consumed by FBS 1 at each boundary caused by K-means clustering of all instances of the monte carlo simulation for every set of nodes
ec1_mean_20nodes_Kmeans = mean(ec1_kmeans(1:montecarlo_i,:));
ec1_mean_30nodes_Kmeans = mean(ec1_kmeans(montecarlo_i+1:2*montecarlo_i,:));
ec1_mean_40nodes_Kmeans = mean(ec1_kmeans(2*montecarlo_i+1:3*montecarlo_i,:));
ec1_mean_50nodes_Kmeans = mean(ec1_kmeans(3*montecarlo_i+1:4*montecarlo_i,:));
ec1_mean_60nodes_Kmeans = mean(ec1_kmeans(4*montecarlo_i+1:5*montecarlo_i,:));
ec1_mean_70nodes_Kmeans = mean(ec1_kmeans(5*montecarlo_i+1:6*montecarlo_i,:));
ec1_mean_80nodes_Kmeans = mean(ec1_kmeans(6*montecarlo_i+1:7*montecarlo_i,:));

%creating matrix containing mean energy consumed by FBS 1 at the optimal boundary caused by K-means clustering for each set of nodes
e1_opt_allnodes_Kmeans = [min(ec1_mean_20nodes_Kmeans) min(ec1_mean_30nodes_Kmeans) min(ec1_mean_40nodes_Kmeans) min(ec1_mean_50nodes_Kmeans) min(ec1_mean_60nodes_Kmeans) min(ec1_mean_70nodes_Kmeans) min(ec1_mean_80nodes_Kmeans)];


%visualization of energy consumed by FBS 2 data---------------------------------------------------------------------------------------------------------------------------------------

%finding the mean energy consumed by FBS 2 at each piecewise boundary of all instances of the monte carlo simulation for every set of nodes
ec2_mean_20nodes = mean(energy_consumption2_pw(1:montecarlo_i,:));
ec2_mean_30nodes = mean(energy_consumption2_pw(montecarlo_i+1:2*montecarlo_i,:));
ec2_mean_40nodes = mean(energy_consumption2_pw(2*montecarlo_i+1:3*montecarlo_i,:));
ec2_mean_50nodes = mean(energy_consumption2_pw(3*montecarlo_i+1:4*montecarlo_i,:));
ec2_mean_60nodes = mean(energy_consumption2_pw(4*montecarlo_i+1:5*montecarlo_i,:));
ec2_mean_70nodes = mean(energy_consumption2_pw(5*montecarlo_i+1:6*montecarlo_i,:));
ec2_mean_80nodes = mean(energy_consumption2_pw(6*montecarlo_i+1:7*montecarlo_i,:));

%creating matrix containing mean energy consumed by FBS 2 at the optimal boundary for each set of nodes 
ec2_opt_allnodes = [min(ec2_mean_20nodes) min(ec2_mean_30nodes) min(ec2_mean_40nodes) min(ec2_mean_50nodes) min(ec2_mean_60nodes) min(ec2_mean_70nodes) min(ec2_mean_80nodes)];

%finding the mean energy consumed by FBS 2 at each boundary caused by K-means clustering of all instances of the monte carlo simulation for every set of nodes
ec2_mean_20nodes_Kmeans = mean(ec2_kmeans(1:montecarlo_i,:));
ec2_mean_30nodes_Kmeans = mean(ec2_kmeans(montecarlo_i+1:2*montecarlo_i,:));
ec2_mean_40nodes_Kmeans = mean(ec2_kmeans(2*montecarlo_i+1:3*montecarlo_i,:));
ec2_mean_50nodes_Kmeans = mean(ec2_kmeans(3*montecarlo_i+1:4*montecarlo_i,:));
ec2_mean_60nodes_Kmeans = mean(ec2_kmeans(4*montecarlo_i+1:5*montecarlo_i,:));
ec2_mean_70nodes_Kmeans = mean(ec2_kmeans(5*montecarlo_i+1:6*montecarlo_i,:));
ec2_mean_80nodes_Kmeans = mean(ec2_kmeans(6*montecarlo_i+1:7*montecarlo_i,:));

%creating matrix containing mean energy consumed by FBS 2 at the optimal boundary caused by K-means clustering for each set of nodes
e2_opt_allnodes_Kmeans = [min(ec2_mean_20nodes_Kmeans) min(ec2_mean_30nodes_Kmeans) min(ec2_mean_40nodes_Kmeans) min(ec2_mean_50nodes_Kmeans) min(ec2_mean_60nodes_Kmeans) min(ec2_mean_70nodes_Kmeans) min(ec2_mean_80nodes_Kmeans)];

%plotting #GNs vs Energy consumed graph of FBS 1 comparing FBS Boundary-search clustering to K-means clustering
figure(4)
plot([20:10:80], ec1_opt_allnodes, 'ks--',[20:10:80], e1_opt_allnodes_Kmeans,'yo-','LineWidth',2) 
legend('FBS Boundary-search','K-means','location','nw')
xlabel('Total Number of GNs')
ylabel('Energy Consumption (Ws)')
title('Effect of Clustering on Energy Consumption of FBS 1')
grid on

%plotting #GNs vs Energy consumed graph of FBS 2 comparing FBS Boundary-search clustering to K-means clustering
figure(5)
plot([20:10:80], ec2_opt_allnodes, 'ks--',[20:10:80], e2_opt_allnodes_Kmeans,'yo-','LineWidth',2) 
legend('FBS Boundary-search','K-means','location','nw')
xlabel('Total Number of GNs')
ylabel('Energy Consumption (Ws)')
title('Effect of Clustering on Energy Consumption of FBS 2')
grid on

%bar chart results for every set of nodes containing mean net distance at each boundary and highlighting the optimal boundary
figure(6)
hold on
for bc1 = 1:length(dnet_mean_20nodes)
    result_20nodes = bar(bc1,dnet_mean_20nodes(bc1));
    if dnet_mean_20nodes(bc1) == min(dnet_mean_20nodes)
        set(result_20nodes,'FaceColor','g');
    else
        set(result_20nodes,'FaceColor','r');
    end
    title('20 Nodes: Optimal Piecewise Boundary')
    xlabel('Boundary x-coordinates of p2,p3 (km)')
    ylabel('Mean Net-Distance (km)')                                                                                                                                     
    set(gca,'xtick',1:25,'XTickLabel',{'1.35,1.35','1.35,1.425','1.35,1.5','1.35,1.575','1.35,1.65','1.425,1.35','1.425,1.425','1.425,1.5','1.425,1.575','1.425,1.65','1.5,1.35','1.5,1.425','1.5,1.5','1.5,1.575','1.5,1.65','1.575,1.35','1.575,1.425','1.575,1.5','1.575,1.575','1.575,1.65','1.65,1.35','1.65,1.425','1.65,1.5','1.65,1.575','1.65,1.65'})
    xtickangle(45)
end
hold off

figure(7)
hold on
for bc2 = 1:length(dnet_mean_30nodes)
    result_30nodes = bar(bc2,dnet_mean_30nodes(bc2));
    if dnet_mean_30nodes(bc2) == min(dnet_mean_30nodes)
        set(result_30nodes,'FaceColor','g');
    else
        set(result_30nodes,'FaceColor','r');
    end
    title('30 Nodes: Optimal Piecewise Boundary')
    xlabel('Boundary x-coordinates of p2,p3 (km)')
    ylabel('Mean Net-Distance (km)')
    set(gca,'xtick',1:25,'XTickLabel',{'1.35,1.35','1.35,1.425','1.35,1.5','1.35,1.575','1.35,1.65','1.425,1.35','1.425,1.425','1.425,1.5','1.425,1.575','1.425,1.65','1.5,1.35','1.5,1.425','1.5,1.5','1.5,1.575','1.5,1.65','1.575,1.35','1.575,1.425','1.575,1.5','1.575,1.575','1.575,1.65','1.65,1.35','1.65,1.425','1.65,1.5','1.65,1.575','1.65,1.65'})
    xtickangle(45)
end
hold off 

figure(8)
hold on
for bc3 = 1:length(dnet_mean_40nodes)
    result_40nodes = bar(bc3,dnet_mean_40nodes(bc3));
    if dnet_mean_40nodes(bc3) == min(dnet_mean_40nodes)
        set(result_40nodes,'FaceColor','g');
    else
        set(result_40nodes,'FaceColor','r');
    end
    title('40 Nodes: Optimal Piecewise Boundary')
    xlabel('Boundary x-coordinates of p2,p3 (km)')
    ylabel('Mean Net-Distance (km)')
    set(gca,'xtick',1:25,'XTickLabel',{'1.35,1.35','1.35,1.425','1.35,1.5','1.35,1.575','1.35,1.65','1.425,1.35','1.425,1.425','1.425,1.5','1.425,1.575','1.425,1.65','1.5,1.35','1.5,1.425','1.5,1.5','1.5,1.575','1.5,1.65','1.575,1.35','1.575,1.425','1.575,1.5','1.575,1.575','1.575,1.65','1.65,1.35','1.65,1.425','1.65,1.5','1.65,1.575','1.65,1.65'})
    xtickangle(45)
end
hold off 

figure(9)
hold on
for bc4 = 1:length(dnet_mean_50nodes)
    result_50nodes = bar(bc4,dnet_mean_50nodes(bc4));
    if dnet_mean_50nodes(bc4) == min(dnet_mean_50nodes)
        set(result_50nodes,'FaceColor','g');
    else
        set(result_50nodes,'FaceColor','r');
    end
    title('50 Nodes: Optimal Piecewise Boundary')
    xlabel('Boundary x-coordinates of p2,p3 (km)')
    ylabel('Mean Net-Distance (km)')
    set(gca,'xtick',1:25,'XTickLabel',{'1.35,1.35','1.35,1.425','1.35,1.5','1.35,1.575','1.35,1.65','1.425,1.35','1.425,1.425','1.425,1.5','1.425,1.575','1.425,1.65','1.5,1.35','1.5,1.425','1.5,1.5','1.5,1.575','1.5,1.65','1.575,1.35','1.575,1.425','1.575,1.5','1.575,1.575','1.575,1.65','1.65,1.35','1.65,1.425','1.65,1.5','1.65,1.575','1.65,1.65'})
    xtickangle(45)
end
hold off

figure(10)
hold on
for bc5 = 1:length(dnet_mean_60nodes)
    result_60nodes = bar(bc5,dnet_mean_60nodes(bc5));
    if dnet_mean_60nodes(bc5) == min(dnet_mean_60nodes)
        set(result_60nodes,'FaceColor','g');
    else
        set(result_60nodes,'FaceColor','r');
    end
    title('60 Nodes: Optimal Piecewise Boundary')
    xlabel('Boundary x-coordinates of p2,p3 (km)')
    ylabel('Mean Net-Distance (km)')
    set(gca,'xtick',1:25,'XTickLabel',{'1.35,1.35','1.35,1.425','1.35,1.5','1.35,1.575','1.35,1.65','1.425,1.35','1.425,1.425','1.425,1.5','1.425,1.575','1.425,1.65','1.5,1.35','1.5,1.425','1.5,1.5','1.5,1.575','1.5,1.65','1.575,1.35','1.575,1.425','1.575,1.5','1.575,1.575','1.575,1.65','1.65,1.35','1.65,1.425','1.65,1.5','1.65,1.575','1.65,1.65'})
    xtickangle(45)
end
hold off

figure(11)
hold on
for bc6 = 1:length(dnet_mean_70nodes)
    result_70nodes = bar(bc6,dnet_mean_70nodes(bc6));
    if dnet_mean_70nodes(bc6) == min(dnet_mean_70nodes)
        set(result_70nodes,'FaceColor','g');
    else
        set(result_70nodes,'FaceColor','r');
    end
    title('70 Nodes: Optimal Piecewise Boundary')
    xlabel('Boundary x-coordinates of p2,p3 (km)')
    ylabel('Mean Net-Distance (km)')
    set(gca,'xtick',1:25,'XTickLabel',{'1.35,1.35','1.35,1.425','1.35,1.5','1.35,1.575','1.35,1.65','1.425,1.35','1.425,1.425','1.425,1.5','1.425,1.575','1.425,1.65','1.5,1.35','1.5,1.425','1.5,1.5','1.5,1.575','1.5,1.65','1.575,1.35','1.575,1.425','1.575,1.5','1.575,1.575','1.575,1.65','1.65,1.35','1.65,1.425','1.65,1.5','1.65,1.575','1.65,1.65'})
    xtickangle(45)
end
hold off

figure(12)
hold on
for bc7 = 1:length(dnet_mean_80nodes)
    result_80nodes = bar(bc7,dnet_mean_80nodes(bc7));
    if dnet_mean_80nodes(bc7) == min(dnet_mean_80nodes)
        set(result_80nodes,'FaceColor','g');
    else
        set(result_80nodes,'FaceColor','r');
    end
    title('80 Nodes: Optimal Piecewise Boundary')
    xlabel('Boundary x-coordinates of p2,p3 (km)')
    ylabel('Mean Net-Distance (km)')
    set(gca,'xtick',1:25,'XTickLabel',{'1.35,1.35','1.35,1.425','1.35,1.5','1.35,1.575','1.35,1.65','1.425,1.35','1.425,1.425','1.425,1.5','1.425,1.575','1.425,1.65','1.5,1.35','1.5,1.425','1.5,1.5','1.5,1.575','1.5,1.65','1.575,1.35','1.575,1.425','1.575,1.5','1.575,1.575','1.575,1.65','1.65,1.35','1.65,1.425','1.65,1.5','1.65,1.575','1.65,1.65'})
    xtickangle(45)
end
hold off

%bar chart results for every set of nodes containing mean net distance at optimal boundary vs at K-means boundary 

dnet_pw_kmeans_20 = [min(dnet_mean_20nodes) dnet_mean_20nodes_Kmeans];
dnet_pw_kmeans_30 = [min(dnet_mean_30nodes) dnet_mean_30nodes_Kmeans];
dnet_pw_kmeans_40 = [min(dnet_mean_40nodes) dnet_mean_40nodes_Kmeans];
dnet_pw_kmeans_50 = [min(dnet_mean_50nodes) dnet_mean_50nodes_Kmeans];
dnet_pw_kmeans_60 = [min(dnet_mean_60nodes) dnet_mean_60nodes_Kmeans];
dnet_pw_kmeans_70 = [min(dnet_mean_70nodes) dnet_mean_70nodes_Kmeans];
dnet_pw_kmeans_80 = [min(dnet_mean_80nodes) dnet_mean_80nodes_Kmeans];

figure(13)
hold on
for bc8 = 1:2
    result_20nodes_Kmeans = bar(bc8, dnet_pw_kmeans_20(bc8));
    if min(dnet_pw_kmeans_20(bc8) == min(dnet_pw_kmeans_20))
        set(result_20nodes_Kmeans,'FaceColor','g');
    else
        set(result_20nodes_Kmeans,'FaceColor','r');
    end
    title('20 Nodes: FBS Boundary-search vs. K-means')
    xlabel('Clustering Algorithm')
    ylabel('Mean Net-Distance (km)')
    set(gca,'xticklabel',{[],[],'FBS Boundary-search',[],'K-means',[],[]})
end
hold off

figure(14)
hold on
for bc9 = 1:2
    result_30nodes_Kmeans = bar(bc9, dnet_pw_kmeans_30(bc9));
    if min(dnet_pw_kmeans_30(bc9) == min(dnet_pw_kmeans_30))
        set(result_30nodes_Kmeans,'FaceColor','g');
    else
        set(result_30nodes_Kmeans,'FaceColor','r');
    end
    title('30 Nodes: FBS Boundary-search vs. K-means')
    xlabel('Clustering Algorithm')
    ylabel('Mean Net-Distance (km)')
    set(gca,'xticklabel',{[],[],'FBS Boundary-search',[],'K-means',[],[]})
end
hold off

figure(15)
hold on
for bc10 = 1:2
    result_40nodes_Kmeans = bar(bc10, dnet_pw_kmeans_40(bc10));
    if min(dnet_pw_kmeans_40(bc10) == min(dnet_pw_kmeans_40))
        set(result_40nodes_Kmeans,'FaceColor','g');
    else
        set(result_40nodes_Kmeans,'FaceColor','r');
    end
    title('40 Nodes: FBS Boundary-search vs. K-means')
    xlabel('Clustering Algorithm')
    ylabel('Mean Net-Distance (km)')
    set(gca,'xticklabel',{[],[],'FBS Boundary-search',[],'K-means',[],[]})
end
hold off

figure(16)
hold on
for bc11 = 1:2
    result_50nodes_Kmeans = bar(bc11, dnet_pw_kmeans_50(bc11));
    if min(dnet_pw_kmeans_50(bc11) == min(dnet_pw_kmeans_50))
        set(result_50nodes_Kmeans,'FaceColor','g');
    else
        set(result_50nodes_Kmeans,'FaceColor','r');  
    end
    title('50 Nodes: FBS Boundary-search vs. K-means')
    xlabel('Clustering Algorithm')
    ylabel('Mean Net-Distance (km)')
    set(gca,'xticklabel',{[],[],'FBS Boundary-search',[],'K-means',[],[]})
end 
hold off

figure(17)
hold on
for bc12 = 1:2
    result_60nodes_Kmeans = bar(bc12, dnet_pw_kmeans_60(bc12));
    if min(dnet_pw_kmeans_60(bc12) == min(dnet_pw_kmeans_60))
        set(result_60nodes_Kmeans,'FaceColor','g');
    else
        set(result_60nodes_Kmeans,'FaceColor','r');   
    end
    title('60 Nodes: FBS Boundary-search vs. K-means')
    xlabel('Clustering Algorithm')
    ylabel('Mean Net-Distance (km)')
    set(gca,'xticklabel',{[],[],'FBS Boundary-search',[],'K-means',[],[]})
end 
hold off

figure(18)
hold on
for bc13 = 1:2
    result_70nodes_Kmeans = bar(bc13, dnet_pw_kmeans_70(bc13));
    if min(dnet_pw_kmeans_70(bc13) == min(dnet_pw_kmeans_70))
        set(result_70nodes_Kmeans,'FaceColor','g');
    else
        set(result_70nodes_Kmeans,'FaceColor','r');
    end
    title('70 Nodes: FBS Boundary-search vs. K-means')
    xlabel('Clustering Algorithm')
    ylabel('Mean Net-Distance (km)')
    set(gca,'xticklabel',{[],[],'FBS Boundary-search',[],'K-means',[],[]})
end
hold off

figure(19)
hold on
for bc14 = 1:2
    result_80nodes_Kmeans = bar(bc14, dnet_pw_kmeans_80(bc14));
    if min(dnet_pw_kmeans_80(bc14) == min(dnet_pw_kmeans_80))
        set(result_80nodes_Kmeans,'FaceColor','g');
    else
        set(result_80nodes_Kmeans,'FaceColor','r');
    end
    title('80 Nodes: FBS Boundary-search vs. K-means')
    xlabel('Clustering Algorithm')
    ylabel('Mean Net-Distance (km)')
    set(gca,'xticklabel',{[],[],'FBS Boundary-search',[],'K-means',[],[]})
end
hold off