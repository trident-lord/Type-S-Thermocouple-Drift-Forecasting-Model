# Type-S Thermocouple Drift Forecasting Model

🚀 AI-assisted, Deep Learning–based Digital Transformation for Drift Quantification in Thermocouples  
📍 Developed as part of the **IASc Summer Research Fellowship 2025** at **CSIR–National Physical Laboratory, New Delhi**

---

## 📌 Overview

This repository contains the complete pipeline for digitising, processing, and forecasting drift in **Type S thermocouples** subjected to repeated high-temperature cycles. The project blends **metrological instrumentation**, **data science**, and **machine learning** to build predictive models that quantify and anticipate sensor degradation over time.

---

## 🎯 Objectives

- Digitise precision thermometry instruments (Fluke 1529E, Fluke 1620A)
- Automate data logging, cold junction compensation, and plateau detection
- Quantify EMF–temperature drift across 10 full thermal cycles (600–1200°C)
- Build predictive ML models (Random Forest, XGBoost, FFNN) for drift forecasting
- Compare NIST vs. certificate-based calibration methods
- Lay the foundation for Digital Twin–based recalibration strategies

---

## 🛠️ Components

### `correction_script.ipynb`
> Applies **Cold Junction Compensation** using real-time reference sensor data. Corrects raw EMF values before temperature conversion.

### `preprocessing.ipynb`
> Detects **thermal plateaus**, interpolates EMF at fixed temperatures, and calculates **temperature drift** across cycles.

### `drift.ipynb`
> Performs **cycle-wise drift analysis**, plots EMF-temperature deviations, and compares NIST vs certificate interpolation methods.

### `hysteresis.ipynb`
> Visualises **hysteresis curves** for heating vs cooling phases in each cycle. Validates thermoelectric stability.

### `cycles.ipynb`
> Generates high-resolution plots for all 10 thermal cycles. Used to compare ramp/dwell consistency and experiment reproducibility.

---

## 🤖 Machine Learning

- Models trained: `Linear Regression`, `Polynomial Regression`, `Random Forest`, `XGBoost`, `Feedforward Neural Network`
- Feature Engineering:
  - Plateau temperature
  - Plateau duration
  - EMF deviations
  - Cycle index
- Best model: `Random Forest` with MAE < **0.0000005°C**
- SHAP analysis used for interpretability

---

## 📊 Results Summary

| Model                | MAE (°C)        | MSE (°C²)       |
|---------------------|-----------------|-----------------|
| Baseline            | 4.79 × 10⁻⁷     | 1.08 × 10⁻¹¹    |
| Linear Regression   | 4.78 × 10⁻⁷     | 1.08 × 10⁻¹¹    |
| Polynomial (deg=3)  | 4.80 × 10⁻⁷     | 1.08 × 10⁻¹¹    |
| **Random Forest**   | **4.60 × 10⁻⁷** | **1.02 × 10⁻¹¹**|
| XGBoost             | 4.78 × 10⁻⁷     | 1.08 × 10⁻¹¹    |
| FFNN                | 1.54 × 10⁻⁴     | 3.51 × 10⁻⁸     |

---

## 🔬 Experimental Setup

- **Thermocouple Type**: Type S (Pt–10%Rh / Pt)
- **Temperature Range**: 600°C → 1200°C
- **Instruments**: Fluke 1529E, Fluke 1620A
- **Cycles Logged**: 10 full thermal cycles (~80+ hours)
- **Logging Interval**: 1 Hz
- **Thermal Dwell Enforcement**: ±0.25°C for ≥ 180 seconds

---

## 📦 Future Work

- Extend dataset to 50+ cycles for Type K/N thermocouples
- Introduce **LSTM-based sequential forecasting** models
- Deploy a real-time Digital Twin for automated recalibration alerts
- Enable cloud dashboard for multi-device monitoring

---

## 📄 Acknowledgements

This project was carried out under the **IASc Summer Research Fellowship 2025**, hosted at the **CSIR–National Physical Laboratory (NPL)** under the guidance of:

- [Dr. Dilip D. Shivagan (Senior Principal Scientist)](https://www.linkedin.com/in/dilip-shivagan-90774352)
- [Mr. Ashish Bhatt(Project Scientist)](https://www.linkedin.com/in/ashish-bhatt-0889a9356)
- Mr. Hansraj Meena
- Dr. Komal Bapna (Senior Scientist)
- Mr. Gaurav Gupta
- [Mr. Anshuman Yadav](https://www.linkedin.com/in/anshuman-yadav-343022333)

---

## 📬 Contact

Feel free to reach out:

**Haran Perumal S L**  
B.Tech CSE, SASTRA Deemed University
📧 [126003096@sastra.ac.in]  
🌐 [LinkedIn](linkedin.com/in/haran-perumal-988513203/)

---

## 📜 License

This project is open-source and available under the MIT License.
