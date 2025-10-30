# 🔧 Astrophotonics & Optical Design Toolkit

*A toolkit for astronomical instrumentation development, featuring spectrograph design, fibre optics analysis, and optical performance modelling. Developed for integral field spectrograph design and optimization.*

---

## 🚀 Quick Start

### MATLAB Examples

```matlab
% Spectrograph parameter optimization
[optimal_params, analysis_data] = spectrograph_parameter_sweep({'Y','J','H'}, ...
    'resolving_power', [5000,5000,5000], 'name', 'MCIFU_5000_950');

% Comprehensive geometric analysis
[performance_metrics, geometric_params] = spectrograph_geometric_analysis(...
    'R_Y', 7880, 's1', 7.3e-6, 'nPix', 2000, 'pix', 18e-6);

% Fibre crosstalk analysis
[crosstalk_results, analysis_data] = fibre_crosstalk_simulator('airy', ...
    'fibre_separation', 25e-6, 'wavelength', 1.55e-6);

% Diffraction limit analysis
[transition_data, performance_metrics] = diffraction_limit_analysis(...
    'grating_density', 650e3, 'beam_size', 14.8e-3, 'f_number', 3.57);
```

---

## 📁 Repository Structure

```
Astrophotonics-Toolkit/
├── 📊 MATLAB/
│   ├── Optical_Geometry/           # Spectrograph layout & analysis
│   │   ├── spectrograph_parameter_sweep.m
│   │   ├── spectrograph_geometric_analysis.m  
│   │   └── diffraction_limit_analysis.m
│   ├── Fibre_Optics/               # Fibre bundle & crosstalk analysis
│   │   └── fibre_crosstalk_simulator.m
│   ├── VPHG_Design/                # Volume Phase Holographic Grating tools
│   │   └── (to be added)
│   └── Data_Processing/            # IFS data handling utilities
│       └── (to be added)
├── 🔍 Zemax_Templates/             # Optical design templates
│   ├── Merit_Functions/            # Optimization operands
│   └── Template_Files/             # Quick-start optical designs
└── 🧪 Examples/
    ├── Spectrograph_Design_Example/
    ├── Crosstalk_Analysis_Example/
    └── Resolving_Power_Tradeoff/
```

---

## 🧰 Tool Categories

### Optical System Analysis

* **spectrograph_parameter_sweep.m** — Multi-band parameter optimization and cross-band consistency
* **spectrograph_geometric_analysis.m** — Comprehensive performance analysis with detector coverage verification
* **diffraction_limit_analysis.m** — Geometric vs diffraction-limited performance transition analysis

### Fibre Optics & IFS

* **fibre_crosstalk_simulator.m** — Multi-model PSF analysis (Airy, Gaussian, dispersed spectra) with pixel integration

### VPH Grating Design *(Planned)*

* **kogelnik_efficiency.m** — VPH grating efficiency calculations
* **vphg_parameter_sweep.m** — *d*, Δ*n*, φ optimisation
* **bragg_condition_solver.m** — Optimal incidence angles

### Zemax Automation *(Planned)*

* **ifs_spectrograph_optimization.MF** — Merit functions for spectroscopic systems
* **multi_config_analysis.zpl** — Multi-grating analysis
* **glass_substitution_tool.zpl** — Material optimisation

---

## 📋 Example Workflows

### 1. Spectrograph Design & Optimization

```
Requirements → Parameter sweep → Geometric design → Diffraction analysis → Performance validation
```

### 2. Fibre System Analysis

```
Bundle geometry → Crosstalk simulation → Detector layout → Performance validation
```

### 3. System Performance Budget

```
Geometric resolving power → Diffraction limit → Transition wavelength → Optimization
```

---

## 🎯 Applications

* Astronomical spectrograph design
* Integral Field Spectroscopy (IFS) systems
* Fibre-fed instrument development
* Optical performance modelling and tolerancing
* Cross-dispersed spectrometer design

---

## 🔬 Theory Background

Tools are based on established physical principles:

* **Geometrical optics** — Spectrograph layout and resolving power
* **Fourier optics** — Diffraction analysis and PSF modelling
* **Statistical optics** — Fibre crosstalk and signal analysis
* **Grating theory** — Dispersion and resolution limits

---

## 📝 License & Citation

This toolkit is available for academic and research use.
If used in publications, please cite:

```bibtex
% Your thesis/dissertation citation
@mastersthesis{author2024spectrograph,
  title={Development of an Integral Field Spectrograph for Exoplanet Science},
  author={Your Name},
  year={2024},
  school={Politecnico di Milano}
}
```

---

## 🤝 Contributing

We welcome contributions and enhancements! Please feel free to:

* Submit issues for bugs or feature requests
* Suggest additional tools or improvements
* Share your own spectrograph design utilities

*Tools developed during MSc thesis work on "Development of an Integral Field Spectrograph for Exoplanet Science" at Politecnico di Milano and INAF - Osservatorio Astronomico di Brera.*

