clear all; close all; clc;

% load the CT scan
filePath = '3 FIRST RUN.nrrd';
ds = nrrd_load(filePath);

% show the full bone volume so we can see what we're working with
bone_viewer(ds);

% let the user pick which bone to analyze
roi = select_roi(ds);

% run fracture detection on the selected region
% 300 is the bone density threshold (based on HU ranges for cortical bone)
results = detect_fractures(ds, roi, 300);

% show all the results: charts, slices, and 3D view
show_results(ds, results);