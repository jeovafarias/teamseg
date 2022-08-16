function seg = graph_cut_segmentation(img, theta, lambda)
% GRAPH_CUT_SEGMENTATION  Segment image with region models theta_1 and
% theta_2 using Graph Cuts
%           
%   PARAMS:
%   - img: Image
%   - theta: Region Distributions
%   - lambda: Balancing paramenter
%
%   RETURNS:
%   - seg: final classification as a matrix of the size of img.
    
    % Compute segmentation graph 
    A = make_graph(img, theta, lambda);
    
    % Compute cut
    N = size(theta, 1);     
    lb = graph_cut(A, N);
    
    % Compute final classification
    seg = reshape(lb, size(img));
end

function A = make_graph(img, theta, lambda)
    np = numel(img); 
    N = size(theta, 1); 
    [r, c, ~] = size(img);
    
    % Compute edges from the images's grid graph
    [nx, ny] = meshgrid(1:r, 1:c);
    ct = nx + (ny - 1) * r;
    le = nx + 1 + (ny - 1) * r;
    so = nx + ny * r;

    edges_i = [reshape(ct(:, 1:end-1), [], 1);
               reshape(ct(1:end-1, :), [], 1)]; 
    edges_j = [reshape(le(:, 1:end-1), [], 1);
               reshape(so(1:end-1, :), [], 1)];
    
    % Compute weights from the images's grid graph
    w_grid = lambda * ones(2 * r * c - r - c, 1);
    
    % Compute edges and weights from the st nodes
    w_st_nodes = zeros(size(img, 2), size(img, 1), N);
    for n = 1:N
        % Concatenate new edges
        edges_i = [edges_i; reshape(ct, [], 1)];
        edges_j = [edges_j; (np + n) * ones(np, 1)];
        
        % Compute st weights
        loglikelihood_n = -log(min(1, theta(n, :) + 1e-8));
        w_st_nodes(:, :, n) = loglikelihood_n(img');
    end
    w_st_nodes = reshape(w_st_nodes, [], 1);
    
    % Concatenate all weights
    weights = [w_grid ; w_st_nodes];
    
    % Create sparse adjacency matrix
    A = sparse(edges_i, edges_j, weights, np + N, np + N);
end
