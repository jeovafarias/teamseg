function [img_seg, et, st] = teamseg(img, radius, N, lambda)
% TEAMSET  Apply TEAMSEG algorithm to estimate appearence models of an
% image and segment it using them
%           
%   PARAMS:
%   - img: Image to be processed
%   - radius: neighbourhood radius
%   - N: Number of regions in the image
%   - lambda: segmentation's balancing paramenter
%
%   RETURNS:
%   - img_seg: matrix of labels corresponding to the segmentation of img  
%   - et, st: estimation and segmentation elapsed times

    etstart = tic;
    % Compute the number of colors in img 
    l = length(unique(img(:)));
    
    % Compute the color moments in img
    [alpha, beta, gamma] = compute_color_moments(img, radius, ceil(l/3));
    
    % Estimate appearence based on the moments
    [theta_hat, ~] = tensor_estimator(alpha, beta, gamma, N);
    et = toc(etstart);
    
    ststart = tic;
    % Compute final segmentation
    img_seg = graph_cut_segmentation(img, theta_hat, lambda);
    st = toc(ststart);
end