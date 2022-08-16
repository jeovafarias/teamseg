addpath(genpath('.'))

%% Algorithm Parameters
radius = 25;
lambda = 1;

%% Experiment Parameters
% Choose between running on the available synthetic or real data
use_synthetic = false;

% To select the sythetic image, set desired GT followed by "rand", "gmm30" 
% or "gmm60" (see paper for details)
sythetic_image_name = "gt4_gmm30.png"; 
% Two options for real images: "polsar.png" or "histology.tif"
real_image_name = "polsar.png"; 

%% Read image and GT, if available, and sets N
if use_synthetic
    img_or = imread(strcat("images/synthetic/", sythetic_image_name));
    img = remove_colors(img_or);
    gt_available = true;
    gt = imread(strcat("images/synthetic/", ...
        extractBefore(sythetic_image_name, "_"), ".png"));
    gt = imresize(gt, size(img), "nearest");
    N = length(unique(gt(:)));
else
    img_or = imread(strcat("images/real/", real_image_name));
    img = image_vectorization(double(img_or), 256);
    gt_available = false;
    
    if strcmp(real_image_name, "polsar.png")
        N = 3;
    else
        N = 4;
    end
end 

%% Applies TEAMSEG for estimation and segmentation
[img_seg, et, st] = teamseg(img, radius, N, lambda);

%% Display results
if gt_available
    n_cols = 3;
    colored_seg = color_segmentation(img_seg, gt);
else
    n_cols = 2;
    colored_seg = color_segmentation(img_seg, img_or);
end 
    
subplot(1, n_cols, 1)
imshow(img_or, [])
title("Original Image")
subplot(1, n_cols, 2)
imshow(colored_seg/256, [])
title(sprintf("Segmentation (ET: %.2f s)", et))

if gt_available
    subplot(1, 3, 3)
    imshow(gt, [])
    title("GT")
end



