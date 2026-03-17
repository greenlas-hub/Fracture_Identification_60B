function results = detect_fractures(ds, roi, threshold)
% detect_fractures  Analyze a selected bone region for fractures.
%   results = detect_fractures(ds, roi, threshold)
%
%   Uses solidity tracking: a healthy bone cross-section is convex
%   (solidity close to 1.0). A fracture creates a dent, dropping solidity.
%   We find the healthy part of the bone, then flag where solidity drops.

    V  = ds.data;
    sp = ds.spacing;

    % crop the volume to just the ROI the user selected
    r1 = roi.rows(1);   r2 = roi.rows(2);
    c1 = roi.cols(1);   c2 = roi.cols(2);
    sub = V(r1:r2, c1:c2, :);
    numSlices = size(sub, 3);
    fprintf('ROI sub-volume: %d x %d x %d\n', size(sub,1), size(sub,2), numSlices);

    % make a bone mask by thresholding
    % anything above the threshold is considered bone
    boneMask = sub > threshold;

    % clean up the mask slice by slice using morphological operations
    % imopen removes small noise specks (erode then dilate)
    % imfill fills holes inside the bone (spongy interior)
    % imclose bridges small gaps in the cortex (dilate then erode)
    seClean = strel('disk', 2); % structuring element for cleaning
    seFill  = strel('disk', 3); % structuring element for closing
    for s = 1:numSlices
        slice = boneMask(:,:,s);
        slice = imopen(slice, seClean);
        slice = imfill(slice, 'holes');
        slice = imclose(slice, seFill);
        boneMask(:,:,s) = slice;
    end

    fprintf('Bone voxels in ROI: %d\n', sum(boneMask(:)));

    % isolate the biggest bone blob so we ignore neighboring bones
    % that might also be in the ROI
    % bwconncomp finds all separate 3D objects in the mask
    cc = bwconncomp(boneMask, 26); % 26 = check all neighbor directions
    biggestSize = 0;
    biggestIdx  = 1;
    for i = 1:cc.NumObjects
        blobSize = length(cc.PixelIdxList{i});
        if blobSize > biggestSize
            biggestSize = blobSize;
            biggestIdx  = i;  % remember which blob was biggest
        end
    end

    % zero out everything except the biggest blob
    cleanMask = false(size(boneMask));
    cleanMask(cc.PixelIdxList{biggestIdx}) = true; % only keep the main bone
    boneMask = cleanMask;
    fprintf('Primary bone: %d voxels\n', sum(boneMask(:)));

    % measure the bone shape on each slice using regionprops
    % solidity = area / convex hull area
    %   healthy bone: ~0.95 (smooth, round)
    %   fractured bone: drops to ~0.5-0.7 (has a notch/dent)
    solidityArr = nan(numSlices, 1);
    areaArr     = nan(numSlices, 1);
    piecesArr   = zeros(numSlices, 1);

    for s = 1:numSlices
        sliceMask = boneMask(:,:,s);

        % skip slices with barely any bone
        if sum(sliceMask(:)) < 30
            continue;
        end

        % get properties of the largest region on this slice
        props = regionprops(sliceMask, 'Area', 'Solidity');
        if isempty(props)
            continue;
        end

        % find the biggest region
        bestArea = 0;
        bestIdx = 1;
        for i = 1:length(props)
            if props(i).Area > bestArea
                bestArea = props(i).Area;
                bestIdx = i;
            end
        end
        % save the solidity and area of the biggest region
        solidityArr(s) = props(bestIdx).Solidity;
        areaArr(s)     = props(bestIdx).Area;

        % count how many separate pieces the bone is in on this slice
        % only count pieces bigger than 30 pixels (ignore tiny specks)
        cc2d = bwconncomp(sliceMask, 8); % 8 = check diagonal neighbors too
        count = 0;
        for i = 1:cc2d.NumObjects
            if length(cc2d.PixelIdxList{i}) >= 30
                count = count + 1;
            end
        end
        piecesArr(s) = count;
    end

    % figure out which slices have bone
    boneSlices = find(~isnan(areaArr) & areaArr > 50);
    fprintf('Bone present on %d slices\n', length(boneSlices));

    % find the healthy shaft: longest streak of slices where
    % solidity > 0.90 and the bone is in one piece
    % a healthy bone should be smooth and connected on each slice
    isHealthy = false(numSlices, 1);
    for s = boneSlices(:)'
        if solidityArr(s) > 0.90 && piecesArr(s) == 1
            isHealthy(s) = true;  % this slice looks healthy
        end
    end

    % walk through and find the longest streak of healthy slices
    % similar to counting a winning streak - keep track of current and best
    bestStreakStart = 0;
    bestStreakLen   = 0;
    currentStart    = 0;
    currentLen      = 0;

    for s = 1:numSlices
        if isHealthy(s)
            if currentLen == 0
                currentStart = s;  % start a new streak
            end
            currentLen = currentLen + 1;  % extend the streak
        else
            % streak broke, check if it was the longest so far
            if currentLen > bestStreakLen
                bestStreakLen   = currentLen;
                bestStreakStart = currentStart;
            end
            currentLen = 0;  % reset for next streak
        end
    end
    % check the last streak in case it goes to the end
    if currentLen > bestStreakLen
        bestStreakLen   = currentLen;
        bestStreakStart = currentStart;
    end

    bestStreakEnd = bestStreakStart + bestStreakLen - 1;
    fprintf('Healthy shaft: slices %d-%d (%d slices)\n', ...
        bestStreakStart, bestStreakEnd, bestStreakLen);

    % baseline solidity = median of the healthy shaft
    shaftSolidity = solidityArr(bestStreakStart:bestStreakEnd);
    baseline = median(shaftSolidity, 'omitnan');
    fprintf('Baseline solidity: %.3f\n', baseline);

    % now find a second healthy streak (on the other side of the fracture)
    % zero out the first streak so we can search for the next longest one
    tempHealthy = isHealthy;
    tempHealthy(bestStreakStart:bestStreakEnd) = false;

    % same streak-finding logic as before
    secondStart = 0;
    secondLen   = 0;
    currentStart = 0;
    currentLen   = 0;

    for s = 1:numSlices
        if tempHealthy(s)
            if currentLen == 0
                currentStart = s;
            end
            currentLen = currentLen + 1;
        else
            if currentLen > secondLen
                secondLen   = currentLen;
                secondStart = currentStart;
            end
            currentLen = 0;
        end
    end
    if currentLen > secondLen
        secondLen   = currentLen;
        secondStart = currentStart;
    end
    secondEnd = secondStart + secondLen - 1;

    % the fracture should be between the two healthy streaks
    % if we found two healthy segments, the gap between them is the fracture
    fractureSlices = [];
    gapPerSlice = [];

    if secondLen >= 5
        fprintf('Second healthy segment: slices %d-%d (%d slices)\n', ...
            secondStart, secondEnd, secondLen);

        % figure out which streak comes first in the volume
        if bestStreakStart < secondStart
            gapStart = bestStreakEnd + 1;
            gapEnd   = secondStart - 1;
        else
            gapStart = secondEnd + 1;
            gapEnd   = bestStreakStart - 1;
        end

        % mark slices in the gap as fracture if solidity dropped enough
        % the threshold is baseline minus 0.10 (10% drop = suspicious)
        for s = gapStart:gapEnd
            if isnan(solidityArr(s))
                continue;  % no bone data on this slice
            end
            if solidityArr(s) < (baseline - 0.10)
                fractureSlices(end+1) = s; %#ok
            end
        end
    else
        % only one healthy segment found
        % look nearby (within 30 slices) for anomalies
        searchStart = max(1, bestStreakStart - 30);
        searchEnd   = min(numSlices, bestStreakEnd + 30);

        for s = searchStart:searchEnd
            if isHealthy(s) || isnan(solidityArr(s))
                continue;  % skip healthy slices and empty slices
            end
            % needs a bigger drop (0.15) and decent bone area to count
            if solidityArr(s) < (baseline - 0.15) && areaArr(s) > 200
                fractureSlices(end+1) = s; %#ok
            end
        end
    end

    % measure the gap between bone fragments on each fracture slice
    % if the bone split into 2+ pieces, measure distance between the two biggest
    gapPerSlice = zeros(1, length(fractureSlices));
    for i = 1:length(fractureSlices)
        s = fractureSlices(i);
        sliceMask = boneMask(:,:,s);
        cc2d = bwconncomp(sliceMask, 8);  % find separate pieces
        props = regionprops(cc2d, 'Centroid', 'Area');

        if length(props) >= 2
            % find the two biggest pieces by area
            areas = zeros(length(props), 1);
            for j = 1:length(props)
                areas(j) = props(j).Area;
            end
            [~, order] = sort(areas, 'descend'); % Orders areas in descending order
            c1 = props(order(1)).Centroid;  % center of piece 1
            c2 = props(order(2)).Centroid;  % center of piece 2

            % distance between centers in mm
            dx = (c1(1) - c2(1)) * sp(2);
            dy = (c1(2) - c2(2)) * sp(1);
            gapPerSlice(i) = sqrt(dx^2 + dy^2);
        end
    end

    % build a fracture mask for the 3D visualization
    % bwconvhull gives us what the bone "should" look like if it were intact
    % the difference between hull and actual bone = the fracture notch
    fracMask = false(size(boneMask));
    for i = 1:length(fractureSlices)
        s = fractureSlices(i);
        sliceMask = boneMask(:,:,s);
        hull = bwconvhull(sliceMask);
        % the gap between hull and actual bone = the fracture notch
        fracMask(:,:,s) = (hull & ~sliceMask) | bwperim(sliceMask);
    end

    % print results
    fprintf('\n--- Results ---\n');
    if isempty(fractureSlices)
        fprintf('No fracture detected.\n');
    else
        fprintf('Fracture detected: slices %d-%d (%d slices)\n', ...
            fractureSlices(1), fractureSlices(end), length(fractureSlices));
        if any(gapPerSlice > 0)
            fprintf('Average gap: %.2f mm\n', mean(gapPerSlice(gapPerSlice > 0)));
        end
        fprintf('Solidity in fracture zone: %.3f (vs baseline %.3f)\n', ...
            median(solidityArr(fractureSlices), 'omitnan'), baseline);
    end

    % pack everything into a struct so show_results can use it
    results.subVolume      = sub;
    results.boneMask       = boneMask;
    results.fracMask       = fracMask;
    results.fractureSlices = fractureSlices;
    results.gapPerSlice    = gapPerSlice;
    results.solidityArr    = solidityArr;
    results.areaArr        = areaArr;
    results.piecesArr      = piecesArr;
    results.baseline       = baseline;
    results.threshold      = threshold;
    results.roi            = roi;
    results.spacing        = sp;
    results.boneSlices     = boneSlices;
    results.isHealthy      = isHealthy;
    results.shaftRange     = [bestStreakStart, bestStreakEnd];
end