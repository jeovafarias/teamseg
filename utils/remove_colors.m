function [img_final, map] = remove_colors(img)
% REMOVE_COLORS Remove from the image's histogram the colors that have do
% not show up in the image.

   map = unique(img);
   img_final = img;
   for i = 1:size(img, 1)
        for j = 1:size(img, 2)
            img_final(i, j) = find(map == img(i, j));
        end
   end
end