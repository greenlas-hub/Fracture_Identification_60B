clear; clc; close all;

file = input("nrrd File to analyze: ","s");

Volume  = boneViewer(file);

function newV = boneViewer(file)
info = nrrdinfo(file);
A = info.SpatialMapping.A;
spacing = [norm(A(1:3,1)),norm(A(1:3,2)),norm(A(1:3,3))];
transform = makehgtform('scale',spacing([2,1,3])); % create the transformation for each of the images to have the same custom spacing

medVol = medicalVolume(file);
intensities = medVol.Voxels;
V = single(intensities);
volshow(V,'Transformation',transform);

%sliceViewer(V);

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
