# 🔧 Astrophotonics & Optical Design Toolkit

*A collection of MATLAB and Zemax utilities for astronomical instrumentation development, featuring VPH grating optimisation, spectrograph design, and optical analysis tools.*

---

## 🚀 Quick Start

### MATLAB Examples

```matlab
% VPH grating efficiency calculation
[efficiency, params] = kogelnik_efficiency(...
    'wavelength', 1.3e-6, 'line_density', 650, 'thickness', 20e-6);

% Spectrograph resolving power analysis  
[R_geo, R_diff] = spectrograph_resolving_power(...
    'f_number', 4.55, 'slit_width', 7.3e-6, 'camera_fl', 0.272);

% Fibre crosstalk simulation
crosstalk = fibre_crosstalk_simulator(...
    'fibre_separation', 25e-6, 'wavelength', 1.55e-6);
```

---

## 🔭 Zemax Integration

* Use provided merit functions for spectroscopic systems
* Multi-configuration editors for grating analysis
* Glass substitution templates for material optimisation

---

## 📁 Repository Structure

```
Astrophotonics-Toolkit/
├── 📊 MATLAB/
│   ├── VPHG_Design/           # Volume Phase Holographic Grating tools
│   ├── Optical_Geometry/      # Spectrograph layout & analysis
│   ├── Fibre_Optics/          # Fibre bundle & crosstalk analysis
│   └── Data_Processing/       # IFS data handling utilities
├── 🔍 Zemax_Templates/
│   ├── Merit_Functions/       # Optimisation operands
│   └── Template_Files/        # Quick-start optical designs
├── 📚 Documentation/
│   ├── Getting_Started.md
│   ├── Theory_Background.md
│   └── API_Reference.md
└── 🧪 Examples/
    ├── VPHG_Optimization_Example/
    ├── Crosstalk_Analysis_Example/
    └── Resolving_Power_Tradeoff/
```

---

## 🛠️ Tool Categories

### VPH Grating Design

* `kogelnik_efficiency.m` — VPH grating efficiency calculations
* `vphg_parameter_sweep.m` — *d*, Δ*n*, φ optimisation
* `multiplexed_efficiency.m` — Stacked grating analysis
* `bragg_condition_solver.m` — Optimal incidence angles

### Optical System Analysis

* `spectrograph_resolving_power.m` — Resolving power vs wavelength calculator
* `diffraction_limit_analysis.m` — Airy disk & sampling analysis
* `anamorphic_magnification.m` — Beam compression calculations
* `spot_size_evolution.m` — PSF wavelength dependence

### Fibre Optics & IFS

* `fibre_crosstalk_simulator.m` — PSF overlap analysis
* `datacube_reconstruction.m` — IFS data processing
* `wavelength_calibration.m` — Spectral calibration tools
* `snr_estimator.m` — Signal-to-noise calculations

### Zemax Automation

* `ifs_spectrograph_optimization.MF` — Merit functions
* `multi_config_analysis.zpl` — Multi-grating analysis
* `glass_substitution_tool.zpl` — Material optimisation
* `footprint_diagram_export.zpl` — Batch detector analysis

---

## 📋 Example Workflows

### 1. VPH Grating Design

```
Parameter sweep → Efficiency optimisation → Multiplexing analysis → Manufacturing specs
```

### 2. Spectrograph Layout

```
Requirements → Geometric design → Diffraction analysis → Sampling verification
```

### 3. Fibre System Analysis

```
Bundle geometry → Crosstalk simulation → Detector layout → Performance validation
```

---

## 🎯 Applications

* Astronomical spectrograph design
* Volume Phase Holographic Grating optimisation
* Integral Field Spectroscopy systems
* Fibre-fed instrument development
* Optical performance modelling

---

## 📖 Documentation

* **Getting Started** — Installation and basic usage
* **Theory Background** — Optical design principles
* **API Reference** — Complete function documentation

---

## 🔬 Theory Background

Tools are based on:

* *Kogelnik’s coupled-wave theory* (VPHG efficiency)
* *Fourier optics* (diffraction analysis)
* *Geometrical optics* (spectrograph design)
* *Statistical optics* (crosstalk modelling)

---

## 📝 License & Citation

This toolkit is available for academic use.
If used in research, please cite:

```bibtex
% Your thesis citation would go here
```

---

## 🤝 Contributing

Feel free to submit issues and enhancement requests for additional tools or improved functionality!

> Tools developed during MSc thesis work on *“Development of an Integral Field Spectrograph for Exoplanet Science”* at **Politecnico di Milano** and **INAF - Osservatorio Astronomico di Brera**.

