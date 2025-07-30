% visualize_uncertainty_results.m
function visualize_uncertainty_results(results)
%
% Creates a multi-panel dashboard to visualize how different input
% parameters affect the speed measurement error.
%

% For easier plotting
results.Abs_Error_ms = abs(results.Error_ms);

figure('Name', 'Uncertainty Analysis Dashboard', 'Position', [50 50 1400 800]);
tiledlayout(2, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

% --- Panel 1: Error vs. Wavelength ---
nexttile;
gscatter(results.Wavelength_m / 1000, results.Error_ms, results.Wave_Speed_ms);
xlabel('Wavelength (km)');
ylabel('Speed Error (m/s)');
title('Error vs. Wavelength');
legend('Location','best');
grid on;
ax = gca; ax.Legend.Title.String = 'Wave Speed (m/s)';

% --- Panel 2: Error vs. Alignment of Wave and Drift ---
nexttile;
% Calculate alignment angle: difference between wave and drift direction
alignment_angle = results.Wave_Direction_deg - results.Drift_Angle_deg;
% Use polar scatter to show directional effects
polarscatter(deg2rad(alignment_angle), results.Abs_Error_ms, 30, results.Drift_Speed_ms, 'filled');
title('Absolute Error vs. Wave-Drift Alignment');
pax = gca;
pax.ThetaZeroLocation = 'top';
pax.ThetaDir = 'clockwise';
pax.ThetaTickLabel = {'0° (Aligned)', '45°', '90° (Perp)', '135°', '180° (Opposed)', '', '', ''};
pax.RLim = [0 max(results.Abs_Error_ms)*1.1];
cb = colorbar;
cb.Label.String = 'Drift Speed (m/s)';


% --- Panel 3: Perceived vs. Ground Truth (Accuracy Plot) ---
nexttile;
gscatter(results.Ground_Truth_Speed_ms, results.Perceived_Speed_ms, results.Drift_Speed_ms);
hold on;
% Plot the ideal y=x line
lims = [min([xlim ylim]) max([xlim ylim])];
plot(lims, lims, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Ideal (y=x)');
xlabel('Ground Truth Speed (m/s)');
ylabel('Perceived Speed (m/s)');
title('Overall Algorithm Accuracy');
legend('Location','best');
grid on;
axis equal;
xlim(lims); ylim(lims);
ax = gca; ax.Legend.Title.String = 'Drift Speed (m/s)';

% --- Panel 4: Error vs. Drift Speed ---
nexttile;
gscatter(results.Drift_Speed_ms, results.Error_ms, results.Wavelength_m / 1000);
xlabel('Drift Speed (m/s)');
ylabel('Speed Error (m/s)');
title('Error vs. Drift Speed');
legend('Location','best');
grid on;
ax = gca; ax.Legend.Title.String = 'Wavelength (km)';

% --- Panel 5: Box Plot of Absolute Error by Drift Speed ---
nexttile;
boxplot(results.Abs_Error_ms, results.Drift_Speed_ms, 'symbol', 'rx');
xlabel('Drift Speed (m/s)');
ylabel('Absolute Speed Error (m/s)');
title('Distribution of Absolute Error');
grid on;

% --- Panel 6: Box Plot of Absolute Error by Wavelength ---
nexttile;
boxplot(results.Abs_Error_ms, results.Wavelength_m / 1000, 'symbol', 'bx');
xlabel('Wavelength (km)');
ylabel('Absolute Speed Error (m/s)');
title('Distribution of Absolute Error');
grid on;

end