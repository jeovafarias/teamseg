function [alpha, beta, gamma] = compute_color_moments(img, radius, max_slices)
% COMPUTE_COLOR_MOMENTS  Estimate image moments (1st, 2nd and 3rd order);
%
%   PARAMS:
%   - img: Image to be processed
%   - radius: neighbourhood radius
%   - max_slices: maximum number of slices in the final gamma tensor
%
%   RETURNS:
%   - alpha: First order moment
%   - beta: Second order moment (pairs)
%   - gamma: Third order moment (triples)
    
    % If max_slices is not provided, use all available slices
    if nargin < 3, max_slices = inf; end
    
    % Get image size and number of colors
    l = length(unique(img));
    [nrow, ncol] = size(img);
    
    % Set minimum color to 1 to complay with matlab's indexing
    if min(min(img)) == 0, img = img + 1; end
    
    % Estimate distribuition of colors
    alpha = hist(img(:), 1:l)' / sum(hist(img(:), 1:l));
    alpha = alpha/sum(alpha);
    
    % Initialize pairwise joint distribuition
    beta = zeros(l, l);

    % Initialize three way joint distribuition
    number_slices =  min(max_slices, l);    
    gamma = zeros(l, l, number_slices);
    
    % Sort colors based on their friquencies on I
    [~, ind] = sort(hist(img(:), 1:l), 2, "descend");
    
    % The following vector is used to speed-up checking if a color is of
    % high enough frequency
    is_high_frequency = zeros(l, 1); 
    for i = 1:l
        is_high_frequency(i) = find(ind == i);
    end
    
    for i = 1:nrow
        for j = 1:ncol
            v1 = img(i,j);
            % pixels on bottom side of square centered at (i,j)
            ip = i+radius;
            if ip <= nrow
                for jp = max(j-radius,1):min(j+radius,ncol)
                     v2 = img(ip,jp);
                     beta(v1,v2) = beta(v1,v2)+1;
                end
            end
            % pixels on right side of square centered at (i,j)
            jp = j+radius;
            if jp <= ncol
                for ip = max(i-radius,1):min(i+radius,nrow)
                    v2 = img(ip,jp);
                    beta(v1,v2) = beta(v1,v2)+1;
                end
            end 
            
            % Pixels on the right and botton side
            ip = min(i+radius,nrow);
            jp = min(j+radius,ncol);
            v2 = img(i, jp);
            v3 = img(ip, j);
            
            % Check if color is of high frequency and if it is, counts its
            % ocurences along with its neighbours
            indx = is_high_frequency(v3);
            if indx <= number_slices
                gamma(v1, v2, indx) = gamma(v1, v2, indx) + 1;
                gamma(v2, v1, indx) = gamma(v2, v1, indx) + 1;
            end
            
            indx = is_high_frequency(v2);
            if indx <= number_slices
                gamma(v3, v1, indx) = gamma(v3, v1, indx) + 1;
                gamma(v1, v3, indx) = gamma(v1, v3, indx) + 1;
            end
            
            indx = is_high_frequency(v1);
            if indx <= number_slices
                gamma(v2, v3, indx) = gamma(v2, v3, indx) + 1;
                gamma(v3, v2, indx) = gamma(v3, v2, indx) + 1;
            end
            
            % Pixels on the left and to side
            ip = max(i-radius,1);
            jp = max(j-radius,1);
            v2 = img(i, jp);
            v3 = img(ip, j);
            
            % Repeat what was done in the previous piece of code for the
            % new set of pixels
            indx = is_high_frequency(v3);
            if indx <= number_slices
                gamma(v1, v2, indx) = gamma(v1, v2, indx) + 1;
                gamma(v2, v1, indx) = gamma(v2, v1, indx) + 1;
            end
                
            indx = is_high_frequency(v2);
            if indx <= number_slices
                gamma(v3, v1, indx) = gamma(v3, v1, indx) + 1;
                gamma(v1, v3, indx) = gamma(v1, v3, indx) + 1;
            end
            
            indx = is_high_frequency(v1);
            if indx <= number_slices
                gamma(v2, v3, indx) = gamma(v2, v3, indx) + 1;
                gamma(v3, v2, indx) = gamma(v3, v2, indx) + 1;
            end
        end
    end
    
    % Normalize beta to be a distribution (gamma does not require that)   
    beta = beta + beta';
    beta = beta/sum(beta(:));
end