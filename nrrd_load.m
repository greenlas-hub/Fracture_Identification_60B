function ds = nrrd_load(filepath)
% nrrd_load  Load an NRRD file and return the volume + spacing info.
%   ds = nrrd_load('scan.nrrd')

    % get metadata from the file header
    info = nrrdinfo(filepath);

    % the affine matrix A has direction + spacing baked into its columns
    % the length (norm) of each column = voxel size for that axis
    A = info.SpatialMapping.A;
    spacing = [norm(A(1:3, 1)), norm(A(1:3, 2)), norm(A(1:3, 3))];

    % get the origin from the last column
    origin = A(1:3, 4)';

    % read the actual voxel data
    V = squeeze(double(nrrdread(filepath)));

    % pack into a struct
    ds.data    = V;
    ds.spacing = spacing;
    ds.size    = size(V);

    fprintf('Loaded: %d x %d x %d   spacing: %.3f x %.3f x %.3f mm\n', ...
        ds.size(1), ds.size(2), ds.size(3), spacing(1), spacing(2), spacing(3));
end