function [performance_metrics, geometric_params] = spectrograph_geometric_analysis(varargin)
% SPECTROGRAPH_GEOMETRIC_ANALYSIS - Comprehensive geometric and diffraction analysis
%
% Performs detailed geometric evaluation of spectrograph performance including:
% - Resolving power vs wavelength analysis
% - Detector coverage and spectral length verification
% - Fibre crosstalk and spot overlap assessment
% - Diffraction-limited spot size evolution
% - Multi-band performance optimization
%
% Inputs (as parameter/value pairs):
%   'R_Y'              - Resolving power band Y (default: 7880)
%   's1'               - Slit width [m] (default: 7.3e-6)
%   'nPix'             - Number of pixels (default: 2000)
%   'pix'              - Pixel pitch [m] (default: 18e-6)
%   'F1'               - Collimator f-number (default: 4.55)
%   'samp'             - Sampling factor (default: 2.5)
%   'fibre_separation' - Fibre core separation [m] (default: 25e-6)
%   'rho_Y'            - Line density band Y [lines/m] (default: 650e3)
%   Wavelength ranges with defaults for Y, J, H bands
%
% Example:
%   [metrics, geometry] = spectrograph_geometric_analysis();
%   [metrics, geometry] = spectrograph_geometric_analysis('R_Y', 9000, 'rho_Y', 700e3);

%% Input parsing and parameter initialization
% =========================================================================
p = inputParser;

% Spectrograph core parameters
addParameter(p, 'R_Y', 7880, @isnumeric);                    % Resolving power band Y
addParameter(p, 's1', 7.3e-6, @isnumeric);                   % Slit width [m]
addParameter(p, 'nPix', 2000, @isnumeric);                   % Number of pixels
addParameter(p, 'pix', 18e-6, @isnumeric);                   % Pixel pitch [m]
addParameter(p, 'F1', 4.55, @isnumeric);                     % Collimator f-number
addParameter(p, 'samp', 2.5, @isnumeric);                    % Sampling factor
addParameter(p, 'fibre_separation', 25e-6, @isnumeric);     % Fibre core separation [m]

% Band wavelength ranges [m]
addParameter(p, 'L_min_Y', 0.960e-6, @isnumeric);
addParameter(p, 'L_max_Y', 1.115e-6, @isnumeric);
addParameter(p, 'L_min_J', 1.150e-6, @isnumeric);
addParameter(p, 'L_max_J', 1.345e-6, @isnumeric);
addParameter(p, 'L_min_H', 1.490e-6, @isnumeric);
addParameter(p, 'L_max_H', 1.780e-6, @isnumeric);

% Grating parameters
addParameter(p, 'rho_Y', 650e3, @isnumeric);                % Line density band Y [lines/m]

parse(p, varargin{:});
params = p.Results;

% Calculate derived parameters
params.L_c_Y = (params.L_min_Y + params.L_max_Y) / 2;       % Central wavelength Y [m]
params.L_c_J = (params.L_min_J + params.L_max_J) / 2;       % Central wavelength J [m]
params.L_c_H = (params.L_min_H + params.L_max_H) / 2;       % Central wavelength H [m]
params.length_det = params.nPix * params.pix;               % Detector length [m]

fprintf('=== Spectrograph Geometric Analysis ===\n');
fprintf('Detector: %d pixels, %.1f mm length\n', params.nPix, params.length_det*1e3);
fprintf('Grating density (Y): %.0f lines/mm\n', params.rho_Y*1e-3);

%% Band Y: Primary Design Analysis
% =========================================================================
fprintf('\n--- Band Y Analysis (Reference Design) ---\n');

% Calculate Littrow angle and geometric parameters
alfa = asin((params.rho_Y * params.L_c_Y) / 2);             % Littrow angle [rad]
alfa_deg = rad2deg(alfa);

% Diffraction-limited beam size calculation
W = params.R_Y * 1.22 / params.rho_Y;                       % Beam width [m]
D1 = W * cos(alfa);                                         % Collimator diameter [m]
f1 = params.F1 * D1;                                        % Collimator focal length [m]

% Camera design based on sampling requirements
s2 = params.samp * params.pix;                              % Required spot size [m]
beta_min_Y = asin(params.rho_Y * params.L_min_Y - sin(alfa)); % Min diffraction angle
D2 = W * cos(beta_min_Y);                                   % Camera diameter [m]
F2 = s2 / (2.44 * params.L_min_Y);                          % Camera f-number
f2 = F2 * D2;                                               % Camera focal length [m]

% Geometric spot size and spectral length
s2_geo = params.s1 * F2 / params.F1;                        % Geometric spot size [m]
s2_geo_px = s2_geo / params.pix;                            % Geometric spot size [px]

beta_max_Y = asin(params.rho_Y * params.L_max_Y - sin(alfa)); % Max diffraction angle
beta_diff_Y = beta_max_Y - beta_min_Y;                      % Angular dispersion [rad]
length_Y = beta_diff_Y * f2;                                % Spectral length [m]
length_Y_px = length_Y / params.pix;                        % Spectral length [px]

fprintf('Littrow angle: %.1f°\n', alfa_deg);
fprintf('Beam geometry: W=%.1fmm, D1=%.1fmm, D2=%.1fmm\n', W*1e3, D1*1e3, D2*1e3);
fprintf('Focal lengths: f1=%.1fmm, f2=%.1fmm\n', f1*1e3, f2*1e3);

%% Band J: Cross-Band Consistency
% =========================================================================
fprintf('\n--- Band J Analysis (Cross-Band Matching) ---\n');

% Maintain geometry, calculate required grating density
rho_J = (2 * sin(alfa)) / params.L_c_J;                     % Required line density [lines/m]
R_J = (W * rho_J) / 1.22;                                   % Achievable resolving power

% Spot size range across band
s2_J_min = F2 * 2.44 * params.L_min_J / params.pix;         % Min spot size [px]
s2_J_max = F2 * 2.44 * params.L_max_J / params.pix;         % Max spot size [px]

% Spectral length calculation
beta_min_J = asin(rho_J * params.L_min_J - sin(alfa));
beta_max_J = asin(rho_J * params.L_max_J - sin(alfa));
diff_J = beta_max_J - beta_min_J;
length_J = diff_J * f2;
length_J_px = length_J / params.pix;

fprintf('Line density: %.0f lines/mm\n', rho_J*1e-3);
fprintf('Resolving power: R=%.0f\n', R_J);
fprintf('Spectral length: %.0f pixels\n', length_J_px);

%% Band H: Cross-Band Consistency
% =========================================================================
fprintf('\n--- Band H Analysis (Cross-Band Matching) ---\n');

rho_H = (2 * sin(alfa)) / params.L_c_H;                     % Required line density [lines/m]
R_H = (W * rho_H) / 1.22;                                   % Achievable resolving power

s2_H_min = F2 * 2.44 * params.L_min_H / params.pix;         % Min spot size [px]
s2_H_max = F2 * 2.44 * params.L_max_H / params.pix;         % Max spot size [px]

beta_min_H = asin(rho_H * params.L_min_H - sin(alfa));
beta_max_H = asin(rho_H * params.L_max_H - sin(alfa));
diff_H = beta_max_H - beta_min_H;
length_H = diff_H * f2;
length_H_px = length_H / params.pix;

fprintf('Line density: %.0f lines/mm\n', rho_H*1e-3);
fprintf('Resolving power: R=%.0f\n', R_H);
fprintf('Spectral length: %.0f pixels\n', length_H_px);

%% Performance Summary and Detector Coverage Verification
% =========================================================================
fprintf('\n=== PERFORMANCE SUMMARY ===\n');
display_performance_summary(params, rho_J, rho_H, R_J, R_H, ...
    length_Y_px, length_J_px, length_H_px, f1, f2, D1, D2);

% Verify detector coverage for all bands
coverage_results = verify_detector_coverage(params, f2, ...
    beta_max_Y, beta_max_J, beta_max_H, ...
    beta_min_Y, beta_min_J, beta_min_H, ...
    params.rho_Y, rho_J, rho_H, alfa);

%% Resolving Power Analysis vs Wavelength
% =========================================================================
fprintf('\n--- Resolving Power Analysis ---\n');
[fig_resolving, resolving_data] = analyze_resolving_power(params, ...
    params.rho_Y, rho_J, rho_H, alfa, f2, s2, coverage_results);

%% Fibre Crosstalk and Spot Overlap Analysis
% =========================================================================
fprintf('\n--- Fibre Crosstalk Analysis ---\n');
[crosstalk_metrics, fig_overlap] = analyze_fibre_crosstalk(params, F2, f2, f1);

%% Spot Size Evolution Analysis
% =========================================================================
fprintf('\n--- Spot Size Analysis ---\n');
fig_spotsize = analyze_spot_size_evolution(params, F2);

%% Maximum Fibre Count Calculation
% =========================================================================
max_fibres = calculate_maximum_fibres(params, f2, f1);

%% Prepare Output Structures
% =========================================================================
performance_metrics = struct(...
    'resolving_power', resolving_data, ...
    'detector_coverage', coverage_results, ...
    'crosstalk_analysis', crosstalk_metrics, ...
    'max_fibres', max_fibres, ...
    'spot_sampling', struct('min_Y', s2_J_min, 'max_H', s2_H_max));

geometric_params = struct(...
    'focal_lengths', struct('f1', f1, 'f2', f2), ...
    'diameters', struct('D1', D1, 'D2', D2), ...
    'grating_densities', struct('rho_Y', params.rho_Y, 'rho_J', rho_J, 'rho_H', rho_H), ...
    'angles', struct('alfa', alfa, 'beta_min_Y', beta_min_Y, 'beta_max_Y', beta_max_Y));

fprintf('\n=== Analysis Complete ===\n');

end

%% Helper Functions
% =========================================================================

function display_performance_summary(params, rho_J, rho_H, R_J, R_H, ...
    lenY_px, lenJ_px, lenH_px, f1, f2, D1, D2)

    fprintf('Grating Densities:\n');
    fprintf('  Band Y: %.0f lines/mm\n', params.rho_Y*1e-3);
    fprintf('  Band J: %.0f lines/mm\n', rho_J*1e-3);
    fprintf('  Band H: %.0f lines/mm\n\n', rho_H*1e-3);
    
    fprintf('Resolving Power:\n');
    fprintf('  Band Y: R=%.0f\n', params.R_Y);
    fprintf('  Band J: R=%.0f\n', R_J);
    fprintf('  Band H: R=%.0f\n\n', R_H);
    
    fprintf('Spectral Length (pixels):\n');
    fprintf('  Band Y: %.0f px\n', lenY_px);
    fprintf('  Band J: %.0f px\n', lenJ_px);
    fprintf('  Band H: %.0f px\n\n', lenH_px);
    
    fprintf('Optical System:\n');
    fprintf('  f1 = %.1f mm, f2 = %.1f mm\n', f1*1e3, f2*1e3);
    fprintf('  D1 = %.1f mm, D2 = %.1f mm\n', D1*1e3, D2*1e3);
end

function coverage = verify_detector_coverage(params, f2, ...
    beta_max_Y, beta_max_J, beta_max_H, ...
    beta_min_Y, beta_min_J, beta_min_H, ...
    rho_Y, rho_J, rho_H, alfa)

    coverage = struct();
    spectral_lengths = [...
        (beta_max_Y - beta_min_Y) * f2, ...  % Band Y length [m]
        (beta_max_J - beta_min_J) * f2, ...  % Band J length [m]
        (beta_max_H - beta_min_H) * f2];     % Band H length [m]
    
    bands = {'Y', 'J', 'H'};
    rhos = [rho_Y, rho_J, rho_H];
    beta_maxs = [beta_max_Y, beta_max_J, beta_max_H];
    
    for i = 1:3
        if spectral_lengths(i) <= params.length_det
            fprintf('✓ Band %s fits on detector (%.0f/%.0f px)\n', ...
                bands{i}, spectral_lengths(i)/params.pix, params.nPix);
            coverage.(['band_', bands{i}]) = 'fits';
        else
            beta_cut = (spectral_lengths(i) - params.length_det) / f2;
            lambda_cut = (sin(beta_maxs(i) - beta_cut) + sin(alfa)) / rhos(i);
            
            fprintf('✗ Band %s exceeds detector (%.0f/%.0f px)\n', ...
                bands{i}, spectral_lengths(i)/params.pix, params.nPix);
            fprintf('  Cutoff wavelength: %.1f nm\n', lambda_cut*1e9);
            coverage.(['band_', bands{i}]) = struct('fits', false, 'cutoff_wavelength', lambda_cut);
        end
    end
end

function [fig, resolving_data] = analyze_resolving_power(params, ...
    rho_Y, rho_J, rho_H, alfa, f2, s2, coverage)

    fig = figure('Name', 'Resolving Power Analysis', 'Position', [100, 100, 900, 600]);
    
    % Band Y resolving power calculation
    lambda_Y = linspace(params.L_min_Y, params.L_max_Y, 500);
    beta_Y = asin(rho_Y * lambda_Y - sin(alfa));
    dL_Y = (cos(beta_Y) * s2) ./ (rho_Y * f2);
    rangeR_Y = lambda_Y ./ dL_Y;
    
    % Band J resolving power calculation
    lambda_J = linspace(params.L_min_J, params.L_max_J, 500);
    beta_J = asin(rho_J * lambda_J - sin(alfa));
    dL_J = (cos(beta_J) * s2) ./ (rho_J * f2);
    rangeR_J = lambda_J ./ dL_J;
    
    % Band H resolving power calculation
    lambda_H = linspace(params.L_min_H, params.L_max_H, 500);
    beta_H = asin(rho_H * lambda_H - sin(alfa));
    dL_H = (cos(beta_H) * s2) ./ (rho_H * f2);
    rangeR_H = lambda_H ./ dL_H;
    
    % Plot resolving power
    hold on;
    plot(lambda_Y * 1e9, rangeR_Y, 'LineWidth', 2, 'Color', [0, 0.4470, 0.7410], 'DisplayName', 'Band Y');
    plot(lambda_J * 1e9, rangeR_J, 'LineWidth', 2, 'Color', [0.4660, 0.6740, 0.1880], 'DisplayName', 'Band J');
    plot(lambda_H * 1e9, rangeR_H, 'LineWidth', 2, 'Color', [0.8500, 0.3250, 0.0980], 'DisplayName', 'Band H');
    
    % Add cutoff markers if needed
    colors = [0.8500, 0.3250, 0.0980];
    if isstruct(coverage.band_Y) && ~coverage.band_Y.fits
        xline(coverage.band_Y.cutoff_wavelength * 1e9, '--', 'Color', colors, 'LineWidth', 1.5, ...
            'DisplayName', 'Band Y Cutoff');
    end
    
    xlabel('Wavelength [nm]', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Resolving Power (R)', 'FontSize', 12, 'FontWeight', 'bold');
    title('Resolving Power vs Wavelength', 'FontSize', 14, 'FontWeight', 'bold');
    legend('show', 'Location', 'southwest');
    grid on;
    set(gca, 'FontSize', 11);
    
    resolving_data = struct(...
        'lambda_Y', lambda_Y, 'R_Y', rangeR_Y, ...
        'lambda_J', lambda_J, 'R_J', rangeR_J, ...
        'lambda_H', lambda_H, 'R_H', rangeR_H);
end

function [crosstalk_metrics, fig] = analyze_fibre_crosstalk(params, F2, f2, f1)
    lambda = linspace(params.L_min_Y, params.L_max_H, 1000);
    
    s2_range_1 = F2 * 1.22 * lambda / params.pix;
    s2_range_2 = -F2 * 1.22 * lambda / params.pix;
    
    M = f2 / f1;
    d_fibres_px = params.fibre_separation / params.pix;
    M_d_fibres = d_fibres_px * M;
    
    fig = figure('Name', 'Fibre Crosstalk Analysis', 'Position', [100, 100, 900, 600]);
    hold on;
    
    x = lambda * 1e9;
    y1 = s2_range_1 - M_d_fibres / 2;
    y2 = s2_range_2 - M_d_fibres / 2;
    fill([x, fliplr(x)], [y2, fliplr(y1)], 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.1, ...
        'DisplayName', 'Overlap Region 1');
    
    y3 = s2_range_1 + M_d_fibres / 2;
    y4 = s2_range_2 + M_d_fibres / 2;
    fill([x, fliplr(x)], [y3, fliplr(y4)], 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.1, ...
        'DisplayName', 'Overlap Region 2');
    
    yline(M_d_fibres/2, 'k--', 'LineWidth', 1, 'DisplayName', 'Fibre Boundary');
    yline(-M_d_fibres/2, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
    
    plot(lambda*1e9, s2_range_1 + M_d_fibres/2, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Upper Envelope');
    plot(lambda*1e9, s2_range_2 + M_d_fibres/2, 'b-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    plot(lambda*1e9, s2_range_1 - M_d_fibres/2, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Lower Envelope');
    plot(lambda*1e9, s2_range_2 - M_d_fibres/2, 'r-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    
    xlim([params.L_min_Y*1e9, params.L_max_H*1e9]);
    xlabel('Wavelength [nm]', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Position [pixels]', 'FontSize', 12, 'FontWeight', 'bold');
    title('Fibre Spot Overlap Analysis', 'FontSize', 14, 'FontWeight', 'bold');
    legend('show', 'Location', 'best');
    grid on;
    
    c1 = s2_range_2 + M_d_fibres/2;
    c2 = s2_range_1 - M_d_fibres/2;
    idx = find(c1 < c2, 1);
    
    if ~isempty(idx) && idx < length(lambda)
        lambda_i = lambda(idx) + (lambda(idx+1)-lambda(idx)) .* ...
            (c2(idx)-c1(idx)) ./ ((c1(idx+1)-c1(idx))-(c2(idx+1)-c2(idx)));
        fprintf('Fibre images overlap starting at %.2f μm\n', lambda_i*1e6);
        overlap_occurrence = true;
    else
        fprintf('Fibre images do not overlap within band\n');
        overlap_occurrence = false;
        lambda_i = NaN;
    end
    
    crosstalk_metrics = struct(...
        'overlap_occurs', overlap_occurrence, ...
        'overlap_wavelength', lambda_i, ...
        'fibre_separation_px', M_d_fibres, ...
        'max_spot_size', max(s2_range_1));
end

function fig = analyze_spot_size_evolution(params, F2)
    lambda_Y = linspace(params.L_min_Y, params.L_max_Y, 300);
    lambda_J = linspace(params.L_min_J, params.L_max_J, 300);
    lambda_H = linspace(params.L_min_H, params.L_max_H, 300);
    
    s2_range_Y = F2 * 2.44 * lambda_Y / params.pix;
    s2_range_J = F2 * 2.44 * lambda_J / params.pix;
    s2_range_H = F2 * 2.44 * lambda_H / params.pix;
    
    fig = figure('Name', 'Spot Size Evolution', 'Position', [100, 100, 900, 600]);
    hold on;
    
    plot(lambda_Y * 1e9, s2_range_Y, 'LineWidth', 2, 'Color', [0, 0.4470, 0.7410], 'DisplayName', 'Band Y');
    plot(lambda_J * 1e9, s2_range_J, 'LineWidth', 2, 'Color', [0.4660, 0.6740, 0.1880], 'DisplayName', 'Band J');
    plot(lambda_H * 1e9, s2_range_H, 'LineWidth', 2, 'Color', [0.8500, 0.3250, 0.0980], 'DisplayName', 'Band H');
    
    yline(params.samp, 'r--', 'LineWidth', 2, 'DisplayName', sprintf('Target Sampling (%.1f px)', params.samp));
    
    xlabel('Wavelength [nm]', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Spot Size [pixels]', 'FontSize', 12, 'FontWeight', 'bold');
    title('Diffraction-Limited Spot Size vs Wavelength', 'FontSize', 14, 'FontWeight', 'bold');
    legend('show', 'Location', 'southeast');
    grid on;
    set(gca, 'FontSize', 11);
end

function max_fibres = calculate_maximum_fibres(params, f2, f1)
    M = f2 / f1;
    d_fibres_px = params.fibre_separation / params.pix;
    M_d_fibres = d_fibres_px * M;
    
    n = 1:20;
    N = 1 + 6 * ((n .* (n - 1)) / 2);
    
    L_lim = 2.2e-3 / params.pix;
    L_tot = (N - 1) .* M_d_fibres;
    
    valid_indices = find(L_tot <= L_lim);
    if ~isempty(valid_indices)
        max_fibres_count = N(valid_indices(end));
        fprintf('Maximum number of fibres: %.0f (hexagonal packing)\n', max_fibres_count);
        max_fibres = max_fibres_count;
    else
        warning('No valid fibre configuration found within slit length constraint');
        max_fibres = 0;
    end
end
