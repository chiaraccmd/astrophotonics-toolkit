function [optimal_params, analysis_data] = spectrograph_parameter_sweep(band_selection, varargin)
% SPECTROGRAPH_PARAMETER_SWEEP - Analyze spectrograph parameters vs grating density
%
% Calculates and visualizes how key spectrograph parameters (focal lengths,
% diameters, beam sizes) vary with grating line density for multiple
% photometric bands. Enables cross-band optimization and parameter selection.
%
% Based on geometric optics principles for astronomical spectrograph design:
% - Grating equation and anamorphic magnification
% - Conservation of etendue
% - Resolving power constraints
%
% Inputs:
%   band_selection - Cell array of bands to analyze: {'Y','J','H'} or subset
%   'name'         - Design name string (default: 'MCIFU_5000_950')
%   'resolving_power' - R values for each band (default: [5000,5000,5000])
%   'lambda_central'  - Central wavelengths [m] (default: [1.0375e-6,1.2475e-6,1.635e-6])
%   'slit_width'      - Entrance slit width [m] (default: 7.3e-6)
%   'spatial_res'     - Spatial resolution per band [m] (default: [36e-6,38e-6,50e-6])
%   'f_number_coll'   - Collimator f-number (default: 4.55)
%   'rho_range'       - Grating density range [lines/m] (default: 400e3 to 1200e3)
%
% Outputs:
%   optimal_params    - Structure with optimized parameters for each band
%   analysis_data     - Full parameter sweep data for further analysis
%
% Example:
%   % Analyze all bands with default parameters
%   [opt_params, data] = spectrograph_parameter_sweep({'Y','J','H'});
%
%   % Custom design with specific resolving power
%   [opt_params, data] = spectrograph_parameter_sweep({'Y','J'}, ...
%       'resolving_power', [6000, 6000], 'name', 'HighRes_Design');

%% Input parsing and parameter setup
% =========================================================================
p = inputParser;
addRequired(p, 'band_selection', @iscell);
addParameter(p, 'name', 'MCIFU_5000_950', @ischar);
addParameter(p, 'resolving_power', [5000, 5000, 5000], @isnumeric);
addParameter(p, 'lambda_central', [1.0375e-6, 1.2475e-6, 1.635e-6], @isnumeric);
addParameter(p, 'slit_width', 7.3e-6, @isnumeric);
addParameter(p, 'spatial_res', [36e-6, 38e-6, 50e-6], @isnumeric);
addParameter(p, 'f_number_coll', 4.55, @isnumeric);
addParameter(p, 'rho_range', linspace(400e3, 1200e3, 5000), @isnumeric);

parse(p, band_selection, varargin{:});
params = p.Results;

% Map band names to indices
band_map = containers.Map({'Y','J','H'}, [1,2,3]);
selected_indices = cellfun(@(x) band_map(x), band_selection);

fprintf('Analyzing spectrograph parameters for bands: %s\n', strjoin(band_selection, ', '));

%% Parameter initialization
% =========================================================================
bands = band_selection;
lambdas = [.960, 1.115; 1.150, 1.345; 1.490, 1.780]; % Band edges [μm]

% Extract parameters for selected bands
R = params.resolving_power(selected_indices);
lambdaC = params.lambda_central(selected_indices);
s1 = params.slit_width;
s2 = params.spatial_res(selected_indices);
F1 = params.f_number_coll;
rhoRange = params.rho_range;

n_bands = length(bands);
n_points = length(rhoRange);

fprintf('Parameter sweep range: %.0f to %.0f lines/mm (%d points)\n', ...
    min(rhoRange)*1e-3, max(rhoRange)*1e-3, n_points);

%% Core parameter calculations
% =========================================================================
% Calculate derived optical parameters based on geometric constraints
ratioF = s2 ./ s1;                    % Focal ratio scaling (anamorphic magnification)
F2 = ratioF .* F1;                    % Camera focal ratio

% Preallocate arrays for computational efficiency
deltaLambda = zeros(1, n_bands);
alfa = zeros(n_bands, n_points);      % Incident angle [rad]
f2 = zeros(n_bands, n_points);        % Camera focal length [m]
D = zeros(n_bands, n_points);         % Beam diameter [m]
f1 = zeros(n_bands, n_points);        % Collimator focal length [m]
W = zeros(n_bands, n_points);         % Projected beam size [m]

for i = 1:n_bands
    % Spectral resolution at central wavelength
    deltaLambda(i) = lambdaC(i) / R(i);
    
    % Grating equation: sin(α) = (ρ * λ)/2 for Littrow configuration
    alfa(i, :) = asin((rhoRange .* lambdaC(i)) / 2);
    
    % Camera focal length from spectral dispersion requirement
    f2(i, :) = (cos(alfa(i, :)) .* s2(i)) ./ (rhoRange .* deltaLambda(i));
    
    % Beam diameter (assuming D1 = D2 for simplicity)
    D(i, :) = f2(i, :) ./ F2(i);
    
    % Collimator focal length
    f1(i, :) = F1 .* D(i, :);
    
    % Projected beam size (anamorphic effect)
    W(i, :) = D(i, :) ./ cos(alfa(i, :));
end

%% Visualization
% =========================================================================
% Create comprehensive plots for each selected band
colors = lines(7); % Color scheme for multiple curves

for i = 1:n_bands
    fig = create_parameter_plot(bands{i}, rhoRange, f1(i,:), f2(i,:), D(i,:), colors);
    
    % Save plots automatically
    save_spectrograph_plot(fig, params.name, bands{i});
end

%% Parameter optimization and cross-band matching
% =========================================================================
% Find optimal grating density for reference band (first in selection)
ref_band_idx = 1;
rho_opt = 650e3; % Initial guess or design requirement

[optimal_params, analysis_data] = optimize_cross_band_parameters(...
    bands, rhoRange, rho_opt, f1, f2, D, alfa, lambdaC, ref_band_idx);

%% Display optimization results
% =========================================================================
display_optimization_results(optimal_params, bands);

% Assign outputs
if nargout > 1
    analysis_data.rho_range = rhoRange;
    analysis_data.f1_matrix = f1;
    analysis_data.f2_matrix = f2;
    analysis_data.D_matrix = D;
    analysis_data.alfa_matrix = alfa;
end

fprintf('Parameter sweep completed successfully.\n');

end

%% Helper functions
% =========================================================================

function fig = create_parameter_plot(band_name, rhoRange, f1_band, f2_band, D_band, colors)
% CREATE_PARAMETER_PLOT - Generate formatted plot for band parameters
    fig = figure('Name', sprintf('Spectrograph Parameters - Band %s', band_name), ...
        'Position', [100, 100, 800, 600]);
    hold on;
    
    % Plot main parameters with consistent styling
    h1 = plot(rhoRange * 1e-3, f2_band * 1e3, 'LineWidth', 2.5, ...
        'Color', colors(1,:), 'DisplayName', 'f_2 (camera)');
    h2 = plot(rhoRange * 1e-3, D_band * 1e3, 'LineWidth', 2.5, ...
        'Color', colors(2,:), 'DisplayName', 'D (beam)');
    h3 = plot(rhoRange * 1e-3, f1_band * 1e3, 'LineWidth', 2.5, ...
        'Color', colors(3,:), 'DisplayName', 'f_1 (collimator)');
    
    xlabel('Grating Density [lines/mm]', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Length [mm]', 'FontSize', 12, 'FontWeight', 'bold');
    title(sprintf('Spectrograph Parameters - Band %s', band_name), ...
        'FontSize', 14, 'FontWeight', 'bold');
    
    legend('show', 'Location', 'best');
    grid on;
    set(gca, 'FontSize', 11);
    
    % Set y-axis to log scale for better dynamic range visualization
    set(gca, 'YScale', 'log');
end

function [optimal_params, analysis_data] = optimize_cross_band_parameters(...
    bands, rhoRange, rho_opt, f1, f2, D, alfa, lambdaC, ref_band_idx)
% OPTIMIZE_CROSS_BAND_PARAMETERS - Find consistent parameters across bands

    optimal_params = struct();
    n_bands = length(bands);
    
    % Reference band optimization
    [~, idx_opt] = min(abs(rhoRange - rho_opt));
    optimal_params.ref_band = bands{ref_band_idx};
    optimal_params.rho_opt(ref_band_idx) = rho_opt;
    optimal_params.f1_opt(ref_band_idx) = f1(ref_band_idx, idx_opt);
    optimal_params.f2_opt(ref_band_idx) = f2(ref_band_idx, idx_opt);
    optimal_params.D_opt(ref_band_idx) = D(ref_band_idx, idx_opt);
    optimal_params.alfa_opt(ref_band_idx) = rad2deg(asin((rho_opt * lambdaC(ref_band_idx)) / 2));
    
    fprintf('\n--- Optimization Results ---\n');
    fprintf('Reference band %s at ρ = %.0f lines/mm:\n', ...
        bands{ref_band_idx}, rho_opt*1e-3);
    fprintf('  f1 = %.2f mm, f2 = %.2f mm, D = %.2f mm, α = %.1f°\n', ...
        optimal_params.f1_opt(ref_band_idx)*1e3, ...
        optimal_params.f2_opt(ref_band_idx)*1e3, ...
        optimal_params.D_opt(ref_band_idx)*1e3, ...
        optimal_params.alfa_opt(ref_band_idx));
    
    % Cross-band matching (maintain consistent camera focal length)
    for i = 1:n_bands
        if i ~= ref_band_idx
            [~, idx_match] = min(abs(f2(i,:) - optimal_params.f2_opt(ref_band_idx)));
            optimal_params.rho_opt(i) = rhoRange(idx_match);
            optimal_params.f1_opt(i) = f1(i, idx_match);
            optimal_params.f2_opt(i) = f2(i, idx_match);
            optimal_params.D_opt(i) = D(i, idx_match);
            optimal_params.alfa_opt(i) = rad2deg(asin((optimal_params.rho_opt(i) * lambdaC(i)) / 2));
            
            fprintf('Band %s matched to f2 = %.2f mm:\n', bands{i}, optimal_params.f2_opt(i)*1e3);
            fprintf('  ρ = %.0f lines/mm, f1 = %.2f mm, D = %.2f mm, α = %.1f°\n', ...
                optimal_params.rho_opt(i)*1e-3, optimal_params.f1_opt(i)*1e3, ...
                optimal_params.D_opt(i)*1e3, optimal_params.alfa_opt(i));
        end
    end
    
    analysis_data.matching_tolerance = std([optimal_params.f2_opt]) / mean([optimal_params.f2_opt]);
end

function save_spectrograph_plot(fig, design_name, band_name)
% SAVE_SPECTROGRAPH_PLOT - Save plots to Results directory
    results_dir = 'Results';
    if ~exist(results_dir, 'dir')
        mkdir(results_dir);
        fprintf('Created directory: %s\n', results_dir);
    end
    
    % Save in multiple formats
    fig_name = sprintf('parameters_%s_%s', design_name, band_name);
    saveas(fig, fullfile(results_dir, [fig_name, '.fig']));
    saveas(fig, fullfile(results_dir, [fig_name, '.png']));
    fprintf('Saved plot: %s\n', fig_name);
end

function display_optimization_results(optimal_params, bands)
% DISPLAY_OPTIMIZATION_RESULTS - Print formatted results to console
    fprintf('\n=== FINAL OPTIMIZED PARAMETERS ===\n');
    fprintf('Band\tρ [l/mm]\tf1 [mm]\tf2 [mm]\tD [mm]\tα [°]\n');
    fprintf('----\t--------\t-------\t-------\t------\t-----\n');
    
    for i = 1:length(bands)
        fprintf('%s\t%.0f\t%.1f\t%.1f\t%.1f\t%.1f\n', ...
            bands{i}, optimal_params.rho_opt(i)*1e-3, ...
            optimal_params.f1_opt(i)*1e3, optimal_params.f2_opt(i)*1e3, ...
            optimal_params.D_opt(i)*1e3, optimal_params.alfa_opt(i));
    end
end
