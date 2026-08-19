% ===================================================================================
% An Effective Method for Digital Rock Reconstruction with Enhanced Pore Connectivity
% Author: Chuanyou Zhou
% Date: 2025.05.13
% Description: Processes raw digital rock data, removes small isolated pores,
%              and connects remaining isolated pore clusters to the main pore.
% ===================================================================================
clear all;
clc;
close all;
dbstop if error

tic;
% ********** Open digital rock file, 0 represents pore, 1 represents matrix, .raw format ********
id = fopen('DigitalRock.raw', 'r', 'b');
data = fread(id, 'uint8');
fclose(id);
matrix = reshape(data, [400, 400, 400]); % Digital rock voxel size 400*400*400

% Connected component analysis: identify and label all connected pore regions
CC = bwconncomp(matrix == 0, 18); % Use 18-connectivity
% Calculate the number of voxels (size) for each connected component
voxelSizes = cellfun(@numel, CC.VoxelIdxList);
minVoxelSize = 26;
smallPores = find(voxelSizes < minVoxelSize);
for i = 1:length(smallPores)
    matrix(CC.VoxelIdxList{smallPores(i)}) = 1; % Set these pores as matrix
end

CC = bwconncomp(matrix == 0, 18);
Sizes = size(matrix);
NewCC = true; % CC will be truly updated each time
% Use while loop to process all isolated pore systems
while length(CC.VoxelIdxList) > 1
    voxelSizes = cellfun(@numel, CC.VoxelIdxList);
    [maxPoreSize, maxPoreIndex] = max(voxelSizes);
    maxPore = CC.VoxelIdxList{maxPoreIndex};

    if NewCC
        isolatedPoreIndex = find((1:length(CC.VoxelIdxList)) ~= maxPoreIndex, 1);
    else
        isolatedPoreIndex = find((1:length(CC.VoxelIdxList)) ~= maxPoreIndex, 2);
        isolatedPoreIndex = isolatedPoreIndex(2);
        NewCC = true;
    end

    if isempty(isolatedPoreIndex)
        break;
    end

    % Get the voxel coordinates of the isolated pore system
    poreVoxels = CC.VoxelIdxList{isolatedPoreIndex};
    [x, y, z] = ind2sub(size(matrix), poreVoxels);

    % Calculate the equivalent center point of the isolated pore system and round
    centerX = ceil(mean(x));
    centerY = ceil(mean(y));
    centerZ = ceil(mean(z));
    % Compare (centerX, centerY, centerZ) with all (x, y, z) combinations
    % Check whether the center point is actually within the isolated pore system
    isCenterInPoreVoxels = any(x == centerX & y == centerY & z == centerZ);
    if ~isCenterInPoreVoxels
        centerX = x(1);
        centerY = y(1);
        centerZ = z(1);
    end

    % Calculate the equivalent radius of the isolated pore system and round
    numVoxels = numel(poreVoxels); % Pore volume
    equivalentRadius = ceil((3 * numVoxels / (4 * pi))^(1 / 3)); % Calculate equivalent radius

    % Find the nearest voxel in the largest connected domain
    nearestPore = find_nearest_pore(Sizes, maxPore, centerX, centerY, centerZ);

    % Connect via a random curve with tortuosity 1.5, line width = 1/2 of equivalent sphere radius
    lineWidth = ceil(equivalentRadius/2);
    if lineWidth < 8
        lineWidth = 8;
    end

    % Perform the connection
    matrix = connect_pores_with_curvature(matrix, ...
        [centerX, centerY, centerZ], nearestPore, lineWidth);

    PreCC = CC;
    % Update the largest pore system after each connection
    CC = bwconncomp(matrix == 0, 18);
    if isequal(PreCC, CC)
        % isolatedPoreIndex = find(1:length(CC.VoxelIdxList) ~= maxPoreIndex, 2);
        % isolatedPoreIndex = isolatedPoreIndex(2);
        NewCC = false;
    end
end
% Save the processed digital rock
saveRawData('DigitalRockNew.raw', matrix);
save('DigitalRockNew.mat', 'matrix');

% Print the porosity change before and after
PoreVolume = sum(data == 0); % Count the number of pores before connection
totalVolume = numel(data);
Porosity = PoreVolume / totalVolume; % Calculate porosity
fprintf('Porosity before connection: %.4f\n', Porosity);

finalPoreVolume = sum(matrix(:) == 0); % Count the number of pores after connection
finalPorosity = finalPoreVolume / totalVolume;
fprintf('Porosity after connection: %.4f\n', finalPorosity);

% Close file
fclose(fileID);
toc;

% Function to find the nearest pore
function nearestPore = find_nearest_pore(matrixSize, maxPore, centerX, centerY, centerZ)
% Get the voxel coordinates of the largest connected domain
[conn_x, conn_y, conn_z] = ind2sub(matrixSize, maxPore);

% Calculate the distance from the equivalent center point to the largest pore voxels and round
distances = ceil(sqrt((conn_x - centerX).^2+(conn_y - centerY).^2 ...
    +(conn_z - centerZ).^2));

% Find the largest pore voxel coordinates corresponding to the minimum distance
[~, idx] = min(distances);
nearestPore = [conn_x(idx), conn_y(idx), conn_z(idx)];
nearestPoreW = nearestPore;
% Get points in the nearby region (centered at the nearest point, search within a certain range)
searchRadius = 15;
nearbyPores = [];

for i = 1:length(conn_x)
    if ceil(sqrt((conn_x(i) - nearestPore(1))^2+(conn_y(i) - nearestPore(2))^2 ...
            +(conn_z(i) - nearestPore(3))^2)) <= searchRadius
        nearbyPores = [nearbyPores; conn_x(i), conn_y(i), conn_z(i)];
    end
end

% Calculate the center point of the nearby points and round
if ~isempty(nearbyPores)
    nearestPore = ceil(mean(nearbyPores, 1));
end
% Check whether nearestPore is in the conn_x, conn_y, conn_z list
inList = any(conn_x == nearestPore(1) & conn_y == nearestPore(2) & ...
    conn_z == nearestPore(3));
if ~inList
    nearestPore = nearestPoreW;
end
end

% Function to connect isolated pores (connection line + line width control)
function matrix = connect_pores_with_curvature(matrix, startPore, endPore, lineWidth)
[x1, y1, z1] = deal(startPore(1), startPore(2), startPore(3));
[x2, y2, z2] = deal(endPore(1), endPore(2), endPore(3));

% Calculate the straight-line distance between two points and round
distance = ceil(sqrt((x2 - x1)^2+(y2 - y1)^2+(z2 - z1)^2));

% Generate the connection line
numPoints = 2 * distance;
ConnPoints = generate_random_curve([x1, y1, z1], [x2, y2, z2], numPoints);

% Generate a connection of width lineWidth along the connection line (using radius would require lineWidth/2)
for i = 1:size(ConnPoints, 1)
    currentPoint = ConnPoints(i, :);
    matrix = create_sphere(matrix, currentPoint(1), currentPoint(2), currentPoint(3), ceil(lineWidth/2));
end
end

% Function to generate the connection line
function GenConnPoints = generate_random_curve(startPoint, endPoint, numPoints)
x = ceil(linspace(startPoint(1), endPoint(1), numPoints));
y = ceil(linspace(startPoint(2), endPoint(2), numPoints));
z = ceil(linspace(startPoint(3), endPoint(3), numPoints));
GenConnPoints = [x', y', z'];
end

% Function to create spherical voxels with a certain width along the path
function matrix = create_sphere(matrix, xc, yc, zc, radius)
[rows, cols, depths] = size(matrix);

for x = max(1, ceil(xc-radius)):min(rows, ceil(xc+radius))
    for y = max(1, ceil(yc-radius)):min(cols, ceil(yc+radius))
        for z = max(1, ceil(zc-radius)):min(depths, ceil(zc+radius))
            if ceil(sqrt((x - xc)^2+(y - yc)^2+(z - zc)^2)) <= radius
                matrix(x, y, z) = 0;
            end
        end
    end
end
end

% Data saving function: save data in .raw format
function saveRawData(filename, data)
fileID = fopen(filename, 'w');
fwrite(fileID, data, 'uint8');
fclose(fileID);
end