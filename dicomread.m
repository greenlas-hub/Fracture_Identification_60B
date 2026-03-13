clear; clc; close all;

info = nrrdinfo("1533508883_5CORONAL-imageTypeDERIVED-PRIMARY-AXIAL-CT_SOM5SPO.nrrd");
A = info.SpatialMapping.A;
spacing = [norm(A(1:3,1)),norm(A(1:3,2)),norm(A(1:3,3))];

medVol = medicalVolume("1533508883_5CORONAL-imageTypeDERIVED-PRIMARY-AXIAL-CT_SOM5SPO.nrrd");
intensities = medVol.Voxels;
V = single(intensities);

[L, C] = imsegkmeans3(V, 5,"NormalizeInput",true);
boneMask = (L == 5);
volshow(boneMask, ...
    'Alphamap','linear', ...
    'Transformation',makehgtform('scale',spacing([2,1,3])));
boneMask = single(boneMask);
bone = medicalVolume(boneMask, medVol.VolumeGeometry);


intensitiesSmooth = imgaussfilt3(bone.Voxels,2);
volshow(intensitiesSmooth,...
    'Alphamap','linear',...
    'Transformation',makehgtform('scale',spacing([2,1,3])));
