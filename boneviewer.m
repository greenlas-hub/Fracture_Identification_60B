clear; clc; close all;

file = input("nrrd File to analyze: ","s");

boneViewer(file);

% Loading the Daga
function boneViewer(file)
info = nrrdinfo(file);  % Getting metadata
A = info.SpatialMapping.A;  % Performs voxel scaling and translation
spacing = [norm(A(1:3,1)),norm(A(1:3,2)),norm(A(1:3,3))];
transform = makehgtform('scale',spacing([2,1,3])); % create the transformation for each of the images to have the same custom spacing

% Load the actual 3D volume into the workspace
medVol = medicalVolume(file);
intensities = medVol.Voxels;  % Obtain raw brightness values
V = single(intensities);
volshow(V,'Transformation',transform);

%sliceViewer(V);

mask = (V >= 300) & (V <= 3000); % based on HU unit for bones in CT scans, adjusted for visibility
newV = zeros(size(V),'like',V); % creates newV with same spacing as V
newV(mask) = 1; % every voxel of bone is white
volshow(newV,... 
        'Alphamap','linear', ...
        'Transformation',transform);

% Refining the Image
intensitiesSmooth = imgaussfilt3(newV,2);
volshow(intensitiesSmooth,...
    'Alphamap','linear',...
    'Transformation',transform);
end
