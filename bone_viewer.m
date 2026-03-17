function bone_viewer(ds)
% bone_viewer  Show the full CT scan as a 3D bone volume.
%   bone_viewer(ds) where ds is a struct from nrrd_load

    V = single(ds.data);
    spacing = ds.spacing;

    % transformation so the volume has correct proportions
    transform = makehgtform('scale', spacing([2 1 3]));

    % only show voxels in the bone density range (300-3000 HU)
    mask = (V >= 300) & (V <= 3000);

    % make a new volume where bone = 1, everything else = 0
    boneVol = zeros(size(V), 'like', V);
    boneVol(mask) = 1;

    % smooth it so the surface looks nicer
    boneSmooth = imgaussfilt3(boneVol, 2);

    % show it
    fig = uifigure('Name', 'Bone Viewer', 'Position', [100 100 900 700]);
    v = viewer3d(fig);
    volshow(boneSmooth, ...
        'Parent', v, ...
        'Alphamap', 'linear', ...
        'Transformation', transform);

    fprintf('Showing bone volume.\n');
end