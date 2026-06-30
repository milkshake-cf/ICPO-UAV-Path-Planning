% Quick script: plot Map 1 (Christmas terrain) and save as PNG
clc; clear; close all;

model = CreateModel();

% Mark start and end points
plot3(model.start(1), model.start(2), model.start(3)+20, ...
    'go', 'MarkerSize', 12, 'MarkerFaceColor', 'g');
plot3(model.end(1), model.end(2), model.end(3)+20, ...
    'ro', 'MarkerSize', 12, 'MarkerFaceColor', 'r');

% Adjust view angle
view(135, 30);

% Title and legend
title('Map 1: Christmas Terrain with 6 Threat Cylinders', 'FontSize', 12);
legend({'Terrain', 'Threat', 'Threat', 'Threat', 'Threat', 'Threat', 'Threat', ...
    'Start', 'End'}, 'Location', 'northeast');

% Save figure
saveas(gcf, 'figures/map1_overview.png');
fprintf('Map saved to figures/map1_overview.png\n');
fprintf('Terrain size: %d x %d\n', model.MAPSIZE_X, model.MAPSIZE_Y);
fprintf('Threats: %d cylinders\n', size(model.threats, 1));
fprintf('Start: (%.0f, %.0f, %.0f)\n', model.start);
fprintf('End: (%.0f, %.0f, %.0f)\n', model.end);
