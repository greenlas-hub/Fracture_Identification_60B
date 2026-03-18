clear; clc; close all force;

file = '3 FIRST RUN.nrrd';

% Load CT volume and compute voxel spacing 
[ds.data ds.spacing] = boneViewer(file);
ds.size = size(ds.data);

% Let the user draw a fracture region on a selected slice 

roi = userDrawFractureRegion(ds.data);
disp(roi)

% run fracture detection on the selected region
% 300 is the bone density threshold (based on HU ranges for cortical bone)
results = detect_fractures(ds.data, ds.spacing, roi, 300);

% show all the results: charts, slices, and 3D view
show_results(ds, results);

function [V, spacing] = boneViewer(file)

    % Read NRRD and extract voxel spacing
    info = nrrdinfo(file);
    A = info.SpatialMapping.A;
    spacing = [norm(A(1:3,1)),norm(A(1:3,2)),norm(A(1:3,3))];

    transform = makehgtform('scale',spacing([2,1,3])); % create the transformation for each of the images to have the same custom spacing

    V = squeeze(double(nrrdread(file))); 

    mask = (V >= 300) & (V <= 3000); % based on HU unit for bones in CT scans, adjusted for visibility 
    newV = zeros(size(V),'like',V); % creates newV with same spacing as V
    newV(mask) = 1; % every voxel of bone is white
    volshow(newV,... 
        'Alphamap','linear', ...
        'Transformation',transform);
end

function roi = userDrawFractureRegion(V)
    
    % Display middle slice for ROI selection 
    figure;
    sliceIndex = round(size(V,3)/2); % Choose central slice 
    slice = V(:,:,sliceIndex);

    % Display slice with contrast normalization 
    imshow(mat2gray(slice, double(prctile(slice(:), [5 95]))));
    title("Draw fracture region");

    % Let user draw a rectangular ROI on the slice 
    h = drawrectangle('Color', 'r', 'LineWidth', 2); 
    wait(h); % Wait until ROI is finalized 

     pos = round(h.Position);
    c1 = max(1, pos(1));
    r1 = max(1, pos(2));
    c2 = min(size(V, 2), pos(1) + pos(3));
    r2 = min(size(V, 1), pos(2) + pos(4));

    close(figure);

    % pack into a struct
    roi.rows   = [r1, r2];
    roi.cols   = [c1, c2];
    roi.slices = [1, size(V, 3)];

    fprintf('ROI selected: rows %d-%d, cols %d-%d (%d x %d pixels)\n', ...
        r1, r2, c1, c2, r2-r1+1, c2-c1+1);

    % Convert the drawn ROI into a binary mask 
    close(gcf); % Close ROI figure 
    
end
