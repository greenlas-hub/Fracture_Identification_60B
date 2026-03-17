function show_results(ds, results)
% show_results  Display fracture detection results.
%   Figure 3: full 3D volume with fracture in red

    sp       = results.spacing;
    boneMask = results.boneMask;
    roi      = results.roi;
    % full 3D volume with fracture highlighted in red ---
    % we render two overlapping volumes in the same viewer:
    %   1) the whole scan in white/gray (bone colormap)
    %   2) just the fractured bone in red
    fprintf('Building 3D view...\n');

    V = double(ds.data);

    % normalize to [0 1] using a bone window (200-1500 HU)
    % anything below 200 = air/soft tissue, above 1500 = very dense
    cLim = [200 1500];
    Vn = (V - cLim(1)) / (cLim(2) - cLim(1));
    Vn = max(0, min(1, Vn));

    % map fracture slices from ROI coords back to full volume coords
    fracSlicesFull = [];
    if ~isempty(results.fractureSlices)
        fracSlicesFull = results.fractureSlices + roi.slices(1) - 1;
    end

    % build a mask of the fractured bone in full-volume coordinates
    % we need to map from ROI coords back to the original volume
    % so only the actual bone gets colored red, not everything in those slices
    roiBoneFull = false(size(V));
    if ~isempty(fracSlicesFull)
        s1 = min(fracSlicesFull);
        s2 = max(fracSlicesFull);
        r1 = roi.rows(1); r2 = roi.rows(2);
        c1 = roi.cols(1); c2 = roi.cols(2);

        % copy the bone mask back into full volume coordinates
        for s = s1:s2
            roiSlice = s - roi.slices(1) + 1;  % convert to ROI index
            if roiSlice >= 1 && roiSlice <= size(results.boneMask, 3)
                roiBoneFull(r1:r2, c1:c2, s) = results.boneMask(:,:,roiSlice);
            end
        end
    end

    % split into two volumes: normal bone and fracture bone
    Vn_bone = Vn;
    Vn_bone(roiBoneFull) = 0;  % remove fracture bone from normal volume

    Vn_frac = zeros(size(Vn));
    Vn_frac(roiBoneFull) = Vn(roiBoneFull);  % only fracture bone here

    % alpha ramp: controls which densities are visible
    % low values = transparent (air, soft tissue)
    % high values = opaque (bone)
    x = linspace(0, 1, 256)';
    amap = zeros(256, 1);
    rampStart = 0.15;
    rampEnd   = 0.55;
    rampMask  = x > rampStart;
    amap(rampMask) = ((x(rampMask) - rampStart) / (rampEnd - rampStart)).^1.8;
    amap = min(amap, 0.85);

    % colormaps: white/gray for normal bone, red for fracture
    boneCmap = bone(256);

    redCmap = zeros(256, 3);
    redCmap(:,1) = linspace(0.2, 1.0, 256)'; % red channel goes from dark to bright
    redCmap(:,2) = linspace(0.0, 0.1, 256)';
    redCmap(:,3) = linspace(0.0, 0.05, 256)';

    % render both volumes in the same 3D viewer
    transform = makehgtform('scale', ds.spacing([2 1 3]));

    fig = uifigure('Name', '3D Fracture View', 'Position', [100 50 900 750]);
    v = viewer3d(fig);

    % render normal bone in white/gray
    volshow(Vn_bone, 'Parent', v, ...
        'RenderingStyle', 'VolumeRendering', ...
        'Colormap', boneCmap, 'Alphamap', amap, ...
        'Transformation', transform);

    % render fracture region bone in red
    if any(Vn_frac(:) > 0)
        volshow(Vn_frac, 'Parent', v, ...
            'RenderingStyle', 'VolumeRendering', ...
            'Colormap', redCmap, 'Alphamap', amap, ...
            'Transformation', transform);
    end

    if ~isempty(fracSlicesFull)
        fprintf('Fracture highlighted: slices %d-%d\n', ...
            fracSlicesFull(1), fracSlicesFull(end));
    end
    fprintf('Done.\n');
end