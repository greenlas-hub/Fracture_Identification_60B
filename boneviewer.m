clear; clc; close all;

file = input("nrrd File to analyze: ","s");

% Load CT volume and compute voxel spacing 
V = boneViewer(file);

% Let the user draw a fracture region on a selected slice 
fractureMask3D = userDrawFractureRegion(V);

% Visualize the 3D fracture mask 
volshow(fractureMask3D, 'Alphamap','linear');

fractureMask3D;

function [V, spacing] = boneViewer(file)

    % Read NRRD and extract voxel spacing
    info = nrrdinfo(file);
    A = info.SpatialMapping.A;
    spacing = [norm(A(1:3,1)),norm(A(1:3,2)),norm(A(1:3,3))];

    transform = makehgtform('scale',spacing([2,1,3])); % create the transformation for each of the images to have the same custom spacing

    V = squeeze(double(nrrdread(file))); 

    %sliceViewer(V);

    mask = (V >= 300) & (V <= 3000); % based on HU unit for bones in CT scans, adjusted for visibility 
    newV = zeros(size(V),'like',V); % creates newV with same spacing as V
    newV(mask) = 1; % every voxel of bone is white
    volshow(newV,... 
        'Alphamap','linear', ...
        'Transformation',transform);

   % Smooth intensity visualization 
    intensitiesSmooth = imgaussfilt3(V,2);
    volshow(intensitiesSmooth,...
         'Alphamap','linear',...
        'Transformation',transform);
end

function fractureMask3D = userDrawFractureRegion(V)

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

    % Convert the drawn ROI into a binary mask 
    mask2D = createMask(h); 

    close(gcf); % Close ROI figure 

    % Insert 2D mask into a 3D volume 
    fractureMask3D = false(size(V)); 
    fractureMask3D(:,:,sliceIndex) = mask2D;
    
    % Show selected 2D mask for confirmation
    figure; 
    imshow(mask2D); 
    title("Selected Fracture Mask (2D)"); 
    
end
