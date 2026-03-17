clear; clc; close all;

file = input("nrrd File to analyze: ","s");

V  = boneViewer(file);

fractureMask = userDrawFractureRegion(V);

volshow(fractureMask, 'Alphamap','linear');

fractureMask3D; 

function [V, spacing] = boneViewer(file)
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

intensitiesSmooth = imgaussfilt3(newV,2);
volshow(intensitiesSmooth,...
    'Alphamap','linear',...
    'Transformation',transform);
end

function fractureMask3D = userDrawFractureRegion(V)

    figure;

    sliceIndex = round(size(V,3)/2);
    slice = V(:,:,sliceIndex); 


    imshow(mat2gray(slice, double(prctile(slice(:), [5 95])))); 
    title("Draw fracture region (circle or freehand)");

    disp("Choose ROI type:");
    disp("1 = Circle");
    disp("2 = Freehand");
    roiType = input("Enter choice: ");

    switch roiType
        case 1
            roi = drawcircle('Color','r','LineWidth',1.5);
        case 2
            roi = drawfreehand('Color','r','LineWidth',1.5);
        otherwise
            error("Invalid selection.");
    end

    % Convert ROI to 2-D mask
    mask2D = roi.createMask();

    % Expand to 3-D mask 
    fractureMask3D = false(size(V));
    fractureMask3D(:,:,sliceIndex) = mask2D;

    figure;
    imshow(mask2D);
    title("Selected Fracture Mask (2D)");
end
