

clc
clear all
clf
rng("default")



% Example positions of 10 agents
positions = rand(5, 2) * 10;
plot(positions(:,1),positions(:,2),'.',MarkerSize=20)
hold on
% Distance matrix calculation
distanceMatrix = squareform(pdist(positions));


r = 100; % Maximum connection distance
Deg_lim = 3; % Limit on degree of agents
P_eig = 1; % Panalty for eigenvalue
P_dist = 1; % Panalty for distance


% Find the optimal spanning tree
[optimalAdjMatrix,Optvalue,eig_list] = findOptimalGraph(positions, distanceMatrix, r, Deg_lim, P_eig, P_dist);

%optimalAdjMatrix=[0 1 1;1 0 1;1 1 0];
graph_plot(positions,optimalAdjMatrix);



G=graph(optimalAdjMatrix);
D=degree(G);
L=laplacian(G);
L=full(L);
eigen_values=eig(L);
lambda2=eigen_values(2)
Optvalue

L2 = diag([5,5,5,5,5]) - ones(5,5);
%eig(L2)
%P_eig*(10 - 5)


%% utility funtions

function graph_plot(positions,adjmatrix)

N=length(adjmatrix(1,:));

hold on

for i=1:N
    for j=1:N
        if adjmatrix(i,j)==1
            plot([positions(i, 1) positions(j, 1)],[positions(i, 2) positions(j, 2)]);
        end
    end
end


end


function [optimalAdjMatrix,Optvalue,eig_list] = findOptimalGraph(positions, distanceMatrix, r, Deg_lim, P_eig, P_dist)

    % Number of agents
    numAgents = size(positions, 1);

    % Create the graph with edges only where distance <= r
    G = graph(distanceMatrix <= r & distanceMatrix > 0);

    % Find all connected subgraphs
    subgraphs = allConnectedSubgraphs(G);

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
        objValue = P_eig*(10 - lambda2) + P_dist*totalDistance;
        
        % Update the best graph if this one has a lower objective value
        if objValue < bestObjValue
            bestObjValue = objValue;
            optimalAdjMatrix = adjMatrix;
        end
    end
    Optvalue = bestObjValue;
    % Display the optimal adjacency matrix
    disp('Optimal Adjacency Matrix:');
    disp(optimalAdjMatrix);

    % Plot the optimal graph
    %figure;
    %Goptimal = graph(optimalAdjMatrix);
    %plot(Goptimal, 'Layout', 'force');
    %title('Optimal Connected Graph');
end

% Function to find all connected subgraphs of a graph G
function subgraphs = allConnectedSubgraphs(G)
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
                subgraphs{end+1} = subG; %#ok<AGROW>
            end
        end
    end
end

% Function to check if a graph is connected
function isConn = isconnected(G)
    % A graph is connected if all nodes are part of a single connected component
    isConn = all(conncomp(G) == 1);
end
