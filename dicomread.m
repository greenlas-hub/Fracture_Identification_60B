clear; clc; close all;

% Loading the Data 
info = nrrdinfo("1533508883_5CORONAL-imageTypeDERIVED-PRIMARY-AXIAL-CT_SOM5SPO.nrrd"); % Getting metadata
A = info.SpatialMapping.A; % Performs voxel scaling & translation  
spacing = [norm(A(1:3,1)),norm(A(1:3,2)),norm(A(1:3,3))];

% Load actual 3D volume into the workspace
medVol = medicalVolume("1533508883_5CORONAL-imageTypeDERIVED-PRIMARY-AXIAL-CT_SOM5SPO.nrrd");
intensities = medVol.Voxels; % Obtain raw brightness values
V = single(intensities); 


[L, C] = imsegkmeans3(V, 5,"NormalizeInput",true);
boneMask = (L == 5);
volshow(boneMask, ...
    'Alphamap','linear', ...
    'Transformation',makehgtform('scale',spacing([2,1,3])));
boneMask = single(boneMask);
bone = medicalVolume(boneMask, medVol.VolumeGeometry);

% Refining the Image 
intensitiesSmooth = imgaussfilt3(bone.Voxels,2);
volshow(intensitiesSmooth,...
    'Alphamap','linear',...
    'Transformation',makehgtform('scale',spacing([2,1,3])));




