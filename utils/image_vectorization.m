function [vectorized_img, l] = image_vectorization(img, min_num_partitions, method)
% BUILD_RANDOM_VECTORIZATION  Apply the color space partitioning
%           
%   PARAMS:
%   - img: Image to be processed
%   - min_num_partitions: Minimum number of partitions in the final partitioning
%   - method: Method to be used ('km' for k-means, 'rand' for random)
%
%   RETURNS:
%   - cl_tree: Final partitioning tree
%   - vectorized_img: New color labels
%   - l: Final number of colors

    if nargin < 1, img = double(imread("cameraman.tif")); end
    if nargin < 2, min_num_partitions = 250; end
    if nargin < 3, method = 'rand'; end
    
    tic
    [r, c, dim] = size(img);
    D = reshape(img, [], dim); l = size(unique(D, 'rows'), 1);
    
    % Old partitioning
%     init_idx = 1:length(D); init_id = 1; init_cl = zeros(r, c);
%     [vectorized_img, l] = recursive_partitioning_hyperplanes(D, min_num_pixels_per_partition, init_idx, init_cl, init_id, method);
    
    % Low variance partitioning
    prt.num_partitions = 1; prt.new_img = zeros(r, c);
    prt.variances = []; prt.pixels_in_partition = []; prt.method = method;
    prt.min_num_partitions = min(min_num_partitions, l);
    
    prt = recursive_partitioning_low_variance(D, prt, 0);
    vectorized_img = prt.new_img;
    clear prt
    
    if nargin == 0
        fprintf('Final number of partitions: %d\n', l - 1);
        fprintf('Processing time: %.4f s\n', toc);
        
        figure
        subplot(1, 2, 1)
        imshow(img/255, []);
        title('Original Image')
        subplot(1, 2, 2)
        imshow(vectorized_img, []);
        title('Vectorized Image')
    end
end

function [img, id] = recursive_partitioning_hyperplanes(D, min_num_pixels_per_partition, idx, img, id, method)   
    D_curr = D(idx, :);
    if (numel(idx) < min_num_pixels_per_partition || norm(D_curr - mean(D_curr)) < 1)
        img(idx) = id;
        id = id + 1;
    else
        switch method
            case 'km'
                [labeling, centers] = kmeans(D_curr, 2);

                idx_1 = idx(find(labeling == 1)');
                idx_2 = idx(find(labeling == 2)');
            case 'rand'
                v = rand(size(D_curr, 2), 1);
                Dv = D_curr*v;
                threshold = mean(Dv);

                idx_1 = idx(find(Dv > threshold)');
                idx_2 = idx(find(Dv <= threshold)');    
            otherwise
                error('Wrong Method')
        end
        
        [img, id] = ...
            recursive_partitioning_hyperplanes(D, min_num_pixels_per_partition, idx_1, img, id, method);
        [img, id] = ...
            recursive_partitioning_hyperplanes(D, min_num_pixels_per_partition, idx_2, img, id, method);
    end
end


function prt_info = recursive_partitioning_low_variance(D, prt_info, id)
    
    function prt_info = add_info(prt_info, id, idx_pixels)
        pixels = D(idx_pixels, :);
        mean_pixel = mean(pixels, 1);
        var = sum(vecnorm(pixels - mean_pixel, 2, 2).^2)/length(idx_pixels);
        
        prt_info.variances(id) = var;
        prt_info.pixels_in_partition{id} = idx_pixels;
        prt_info.new_img(idx_pixels) = id;
    end

    if (prt_info.num_partitions < prt_info.min_num_partitions)
        if id == 0
            id = 1; curr_pixels_indexes = 1:size(D, 1);
        else
            curr_pixels_indexes = prt_info.pixels_in_partition{id};
        end
        
        prt_info.num_partitions = prt_info.num_partitions + 1;
        id_child_1 = id; id_child_2 = prt_info.num_partitions;
        
        pixels_curr = D(curr_pixels_indexes, :);
        v = rand(size(pixels_curr, 2), 1); 
        Dv = pixels_curr * v ; threshold = mean(Dv);
        idx_child_1 = curr_pixels_indexes(find(Dv > threshold)');
        idx_child_2 = curr_pixels_indexes(find(Dv <= threshold)');
        
        prt_info = add_info(prt_info, id_child_1, idx_child_1);
        prt_info = add_info(prt_info, id_child_2, idx_child_2);
        [~, id_next] = max(prt_info.variances);
        prt_info = recursive_partitioning_low_variance(D, prt_info, id_next);
    end
end

