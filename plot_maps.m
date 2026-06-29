%_________________________________________________________________________%
%  Plot all 8 benchmark maps as 2x4 subplot grid                           %
%_________________________________________________________________________%

clc; clear; close all;

figure('Position', [100 100 1600 800]);

for m = 1:8
    subplot(2, 4, m);
    try
        model = CreateModel_map(m);
        % Simple 3D plot: terrain + threats + start/end
        surf(model.X, model.Y, model.H, 'EdgeColor', 'none', 'FaceAlpha', 0.7);
        colormap('terrain'); hold on;
        
        % Threats as red cylinders
        for j = 1:size(model.threats, 1)
            [cx, cy, cz] = cylinder(model.threats(j,4), 20);
            cx = cx + model.threats(j,1);
            cy = cy + model.threats(j,2);
            cz = cz * model.threats(j,3) + model.H(round(model.threats(j,2)), round(model.threats(j,1)));
            surf(cx, cy, cz, 'FaceColor', 'r', 'FaceAlpha', 0.4, 'EdgeColor', 'none');
        end
        
        % Start (green) and End (blue) markers
        plot3(model.start(1), model.start(2), model.start(3), 'go', 'MarkerSize', 12, 'MarkerFaceColor', 'g');
        plot3(model.end(1), model.end(2), model.end(3), 'bo', 'MarkerSize', 12, 'MarkerFaceColor', 'b');
        
        title(sprintf('Map %d', m), 'FontSize', 12);
        xlabel('X'); ylabel('Y'); zlabel('Z');
        view(45, 30); axis equal tight;
        
    catch ME
        title(sprintf('Map %d: ERROR', m), 'Color', 'r');
    end
end

sgtitle('8 Benchmark Scenarios for UAV Path Planning', 'FontSize', 14, 'FontWeight', 'bold');
saveas(gcf, 'figures/all_maps.png');
fprintf('Saved to figures/all_maps.png\n');
