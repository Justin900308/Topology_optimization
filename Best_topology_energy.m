%% This script will generate the graph that minimizes the total energy comsumtion
% e is proportional to total distance/lambda2

clc
clear all
clf
rng("default")

N=8;
subplot(2,1,1)
% Example positions of 10 agents
positions = rand(N, 2) * 10;
plot(positions(:,1),positions(:,2),'.',MarkerSize=20)
hold on
% Distance matrix calculation
distanceMatrix = squareform(pdist(positions));


r = 8; % Maximum connection distance
Deg_lim = 5; % Limit on degree of agents



% Find the optimal spanning tree
[optimalAdjMatrix,Optvalue,eig_list] = findOptimalGraph_1(positions, distanceMatrix, r, Deg_lim);

%optimalAdjMatrix=[0 1 1;1 0 1;1 1 0];
graph_plot(positions,optimalAdjMatrix);



G=graph(optimalAdjMatrix);
D=degree(G);
inciden1=incidence(G);
[n,m] = size(inciden1);
w1=ones(1,m)/sum(ones(1,m));
L1=inciden1*diag(w1)*inciden1';
L1=full(L1);
eigen_values=eig(L1);
lambda2=eigen_values(2)




subplot(2,1,2)
hold on
[w2,inciden2] = findOptimalGraph_2(optimalAdjMatrix);
G = graph(optimalAdjMatrix);
G.Edges.Weight = w2;
optimalAdjMatrix2 = adjacency(G,w2);
optimalAdjMatrix2 = full(optimalAdjMatrix2);
graph_plot(positions,optimalAdjMatrix2);
L2=inciden2 * diag(w2) * inciden2';
plot(positions(:,1),positions(:,2),'.',MarkerSize=20)
eigen_values2=eig(L2);
lambda2=eigen_values2(2)
%% 


%% utility funtions


function  [w,inciden]=findOptimalGraph_2(optimalAdjMatrix)

G = graph(optimalAdjMatrix);
A = incidence(G);
[n,m] = size(A);
w=ones(m,1);
%L = A * diag(w) * A';
sum1=0;
%eig1=eig(L);


I = eye(n,n);
J = I - (1/n) * ones(n,n);
cvx_begin sdp 
cvx_precision best 
    variable w(m,1)   % edge weights
    variable r        % epigraph variable
    variable b        % epigraph variable
    variable L(n,n) symmetric
    minimize(10-r)
    %maximize(1*r)
    subject to
        L == A * diag(w) * A';
        r * I <= L+b*ones(n,n) ;
        r>=0;
        for i=1:4
            w(i)>=0;
            w(i)<=1;
        end
        ones(1,m)*w<=1;
cvx_end
inciden=A;


end






function graph_plot(positions,adjmatrix)

N=length(adjmatrix(1,:));
total_weight=sum(sum(adjmatrix));
adjmatrix_normal=adjmatrix/total_weight;
hold on

for i=1:N
    for j=1:N
        if adjmatrix(i,j) > 0
            plot([positions(i, 1) positions(j, 1)],[positions(i, 2) positions(j, 2)], ...
                LineWidth=adjmatrix_normal(i,j)^4*100000);
        end
    end
end


end


function [optimalAdjMatrix,Optvalue,eig_list] = findOptimalGraph_1(positions, distanceMatrix, r, Deg_lim)

    % Number of agents
    numAgents = size(positions, 1);

    % Create the graph with edges only where distance <= r
    G = graph(distanceMatrix <= r & distanceMatrix > 0);

    % Find all connected subgraphs
    subgraphs = allConnectedSubgraphs(G,Deg_lim);

    % Initialize variables to track the best graph
    bestObjValue = inf;
    optimalAdjMatrix = zeros(numAgents);
    eig_list = zeros(length(subgraphs),1);

    % Loop through each connected subgraph to evaluate the objective function
    for k = 1:length(subgraphs)
        % Get the adjacency matrix of the current subgraph
        adjMatrix = adjacency(subgraphs{k});
        adjMatrix = full(adjMatrix);
        % Check if the degree constraint is satisfied (max degree of any node is <= 2)
        degrees = sum(adjMatrix, 2);
        if any(degrees > Deg_lim) || length(degrees) < numAgents
            continue; % Skip subgraphs that violate the degree constraint
            eig_list(k) = 0;
        end
        
        % Compute the Laplacian matrix L = D - A (D is the degree matrix)
        degreeMatrix = diag(sum(adjMatrix, 2));
        laplacianMatrix = degreeMatrix - adjMatrix;
        
        % Compute the eigenvalues of the Laplacian matrix
        eigValues = eig(laplacianMatrix);
        lambda2 = eigValues(2); % Second smallest eigenvalue (algebraic connectivity)
        eig_list(k) = lambda2;
        % Compute the total distance for the current subgraph
        totalDistance = sum(sum(adjMatrix .* distanceMatrix)) / 2; % Total distance
        
        % Objective function: (10 - lambda2) + totalDistance
        % objValue = P_eig*(10 - lambda2) + P_dist*totalDistance;
        objValue = totalDistance / lambda2;
        % Update the best graph if this one has a lower objective value
        if objValue < bestObjValue
            bestObjValue = objValue;
            optimalAdjMatrix = adjMatrix;
        end
    end
    Optvalue = bestObjValue;
    % Display the optimal adjacency matrix
    %disp('Optimal Adjacency Matrix:');
    %disp(optimalAdjMatrix);

    % Plot the optimal graph
    %figure;
    %Goptimal = graph(optimalAdjMatrix);
    %plot(Goptimal, 'Layout', 'force');
    %title('Optimal Connected Graph');
end

% Function to find all connected subgraphs of a graph G
function subgraphs = allConnectedSubgraphs(G, Deg_lim)
    % Get the number of nodes
    n = numnodes(G);
    
    % Initialize an empty list to store all connected subgraphs
    subgraphs = {};
    
    % Extract the list of edges (as pairs of node indices) from the graph G
    edges = G.Edges.EndNodes; % This returns the edge list
    
    % Get the total number of edges in the graph
    numEdges = size(edges, 1);
    
    % Loop over all combinations of edges (from n-1 to numEdges)
    for numSelectedEdges = (n-1):numEdges
        combs = nchoosek(1:numEdges, numSelectedEdges);
        
        % Check each combination of edges to see if it forms a valid connected subgraph
        for i = 1:size(combs, 1)
            % Get the selected edges for this combination
            selectedEdges = edges(combs(i, :), :);
            
            % Create a subgraph using the selected edges
            subG = graph(selectedEdges(:, 1), selectedEdges(:, 2));
            

            % Check if this subgraph is connected
            if isconnected(subG)
                subgraphs{end+1} = subG; 
            end
        end
    end
end

% Function to check if a graph is connected
function isConn = isconnected(G)
    % A graph is connected if all nodes are part of a single connected component
    isConn = all(conncomp(G) == 1);
end
