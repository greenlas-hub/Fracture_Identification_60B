function roi = select_roi(ds)
% select_roi  Let the user draw a box around the bone they want to analyze.
%   roi = select_roi(ds)
%   Shows a mid-slice, user draws a rectangle, double-clicks to confirm.

    V = ds.data;

    % show the middle slice of the volume
    mid = round(size(V, 3) / 2);
    slice = V(:,:,mid);

    fig = figure('Name', 'Draw a box around the bone of interest', ...
        'NumberTitle', 'off', 'Position', [100 100 800 600]);
    imshow(mat2gray(slice, double(prctile(slice(:), [5 95]))));
    title(sprintf('Slice %d - draw a rectangle, then double-click inside it', mid), ...
        'FontSize', 12);

    % let user draw a rectangle
    h = drawrectangle('Color', 'r', 'LineWidth', 2);

    % wait for double-click to confirm
    wait(h);

    % get the position [x y width height] = [col row width height]
    pos = round(h.Position);
    c1 = max(1, pos(1));
    r1 = max(1, pos(2));
    c2 = min(size(V, 2), pos(1) + pos(3));
    r2 = min(size(V, 1), pos(2) + pos(4));

    close(fig);

    % pack into a struct
    roi.rows   = [r1, r2];
    roi.cols   = [c1, c2];
    roi.slices = [1, size(V, 3)];

    fprintf('ROI selected: rows %d-%d, cols %d-%d (%d x %d pixels)\n', ...
        r1, r2, c1, c2, r2-r1+1, c2-c1+1);
end