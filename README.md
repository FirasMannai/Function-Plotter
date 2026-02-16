
# 📈 Function Plotter 🚀

A powerful, multi-mode mathematical function plotter built with **Python**, **Tkinter**, and **Matplotlib**. Whether you need to visualize simple 2D curves, complex 3D surfaces, or analyze control systems with Bode plots, this tool has you covered!

---

## ✨ Features

### 🖼️ Multi-Mode Visualization
* **2D Plotting**: Graph multiple functions simultaneously with customizable colors, line styles, and widths.
* **🌀 Parametric Plots**: Visualize complex curves defined by $x(t)$ and $y(t)$ equations.
* **📶 Bode Diagrams**: Analyze frequency responses (Magnitude & Phase) simply by providing transfer function coefficients.
* **🏔️ 3D Surfaces & Contours**: Explore functions of two variables $f(x, y)$ with interactive 3D rotations and side-by-side contour maps.

---

### 🔬 Advanced Analysis
* **Automated Detection**: Automatically find and mark **Roots**, **Extrema**, and **Inflection Points** on your graphs.
* **Calculus Suite**: Calculate and overlay **1st and 2nd derivatives** on the fly to analyze function behavior visually.

---

### 📊 Data & Integration
* **Data Overlay**: Import `.txt` or `.csv` data files to compare real-world experimental data against mathematical models.
* **💾 Config Management**: Save your entire session settings to a file and reload them later to pick up exactly where you left off.

---

### 🛠️ Flexible Usage
* **Interactive GUI**: Use the intuitive dashboard for real-time exploration and tweaking.
* **Headless Mode**: Generate plots instantly from the command line using configuration files—perfect for automation.

---

## 🚀 Getting Started

### 📋 Prerequisites

Ensure you have Python installed, then install the required dependencies:

```bash
pip install numpy matplotlib

```

### 🏃 Running the App

1. **GUI Mode** (Default):
```bash
python app.py

```


This launches the interactive window where you can tweak all settings.


2. **Headless Mode**:
```bash
python app.py config_bode.txt

```


Directly generates a plot window based on the provided configuration file.



---

## 🛠️ Configuration File Format

The application uses a simple `key: value` format for configuration files.

**Example (`config_bode.txt`):**

```text
# First-order low-pass filter analysis
bode_num: [1]
bode_den: [1, 2, 1]
wmin: 0.1
wmax: 300
title: Bode Diagram – Low-pass Filter

```

---

## 📂 Project Structure

* `app.py`: The main entry point of the application.

* `gui.py`: Handles the Tkinter interface and user interactions.

* `plotter.py`: The engine that performs math evaluations and Matplotlib rendering.

* `config_parser.py`: Logic for reading and converting configuration text files.

* `widgets.py`: Custom UI components like scrollable frames.

---

## 🧩 Supported Math Functions

You can use standard math notation in the input fields, including:

* `sin`, `cos`, `tan`, `exp`, `log`, `log10`, `sqrt`, `abs`
* Constants: `pi`, `e`
* Power operator: `^` (e.g., `x^2`)

---
