function colored_seg = color_segmentation(img_seg, img)
colored_seg = zeros(size(img));
N = length(unique(img_seg(:)));
for c = 1:size(img, 3)
    img_c = img(:, :, c);
    img_seg_c = zeros(size(img, [1, 2]));
    for k =1:N
        img_seg_c(img_seg == k) = mean(img_c(img_seg == k), 'all');
    end
    colored_seg(:, :, c) = img_seg_c;
end
end