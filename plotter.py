
# plotter.py
# All plotting logic (2D, Parametric, Bode, 3D).

import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 (enable 3D backend)

class FunctionPlotter:
    """Implements all plotting modes: 2D, Parametric, Bode, 3D.

    This class contains the core logic for evaluating mathematical expressions
    and generating plots using Matplotlib. It is designed to be used either
    by a GUI or in a headless mode.
    """
    def __init__(self):
        """Initializes the plotter, setting figure and axes to None."""
        self.fig = None
        self.ax = None

    # ----- Public API -----
    def plot_on_canvas(self, figure, ax, config):
        """Draws a plot on a given Matplotlib Figure and Axes.

        This is the main plotting entry point used by the GUI. It orchestrates
        the plotting based on the provided configuration dictionary, handling
        different plot modes (2D, parametric, Bode, 3D).

        Args:
            figure (matplotlib.figure.Figure): The figure to draw on.
            ax (matplotlib.axes.Axes): The axes to draw on. The method may
                replace this if a different plot layout is needed.
            config (dict): A dictionary of configuration parameters.

        Returns:
            matplotlib.axes.Axes: The active axes object after plotting. This
            might be a new axes object if the layout changed.
        """
        # Ensure we have a valid Axes
        if ax is None or ax not in figure.axes:
            ax = figure.axes[0] if figure.axes else figure.add_subplot(111)

        # Special types first (they clear the figure and set own layout)
        if self._maybe_plot_bode(figure, config):
            self.ax = figure.axes[0]
            return self.ax
        if self._maybe_plot_surface_or_contour(figure, config):
            self.ax = figure.axes[0]
            return self.ax

        # Standard single-axes layout
        if len(figure.axes) != 1 or ax not in figure.axes:
            figure.clear()
            ax = figure.add_subplot(111)

        self.ax = ax
        self.ax.clear()

        # 1) x sampling
        xstart = float(config.get('xstart', 0))
        xend   = float(config.get('xend', 5))
        xstep  = float(config.get('xstep', 0.1))
        if xstep <= 0:
            raise ValueError("xstep must be positive.")
        x = np.arange(xstart, xend + xstep, xstep, dtype=float)

        # Parametric?
        if config.get('px') is not None and config.get('py') is not None:
            self._plot_parametric(config)
            self._configure_plot(config)
            return self.ax

        # 2) functions list (string or list)
        functions = config.get('functions', [])
        if isinstance(functions, str):
            functions = [functions]

        # 3) styles (single value or lists accepted; comma-separated allowed)
        colors = config.get('colors') or [config.get('color', 'blue')]
        if isinstance(colors, str):
            colors = [c.strip() for c in colors.split(',')] if ',' in colors else [colors]
        linewidths = config.get('linewidths') or [config.get('linewidth', 2)]
        if isinstance(linewidths, (int, float)):
            linewidths = [linewidths]
        elif isinstance(linewidths, str):
            linewidths = [float(v.strip()) for v in linewidths.split(',')]
        linestyles = config.get('linestyles') or ['-']
        if isinstance(linestyles, str):
            linestyles = [s.strip() for s in linestyles.split(',')] if ',' in linestyles else [linestyles]

        # 4) plot curves + optional features
        for i, func_str in enumerate(functions):
            y = self._evaluate_function(func_str, x)
            c  = colors[i % len(colors)]
            lw = linewidths[i % len(linewidths)]
            ls = linestyles[i % len(linestyles)]
            label = "f(x)" if len(functions) == 1 else f"f{i+1}(x)"
            self.ax.plot(x, y, color=c, linewidth=lw, linestyle=ls, label=label)

            # derivatives
            if config.get('derivative1'):
                dy = self._central_diff(y, xstep)
                self.ax.plot(x, dy, linewidth=max(1.0, lw*0.9), linestyle='--', label=f"{label}'")
            if config.get('derivative2'):
                d2y = self._central_diff(self._central_diff(y, xstep), xstep)
                self.ax.plot(x, d2y, linewidth=max(1.0, lw*0.9), linestyle=':', label=f"{label}''")

            # features
            if any([config.get('find_roots'), config.get('find_extrema'), config.get('find_inflections')]):
                self._find_and_mark_features(x, y, xstep, config)

        # 5) optional data overlay
        if config.get('data_file'):
            self._plot_data_file(config)

        # 6) cosmetics
        self._configure_plot(config)
        return self.ax

    def plot_from_config(self, config):
        """Creates a new plot in a new window from a configuration dictionary.

        This method is used for the headless mode of the application.

        Args:
            config (dict): A dictionary of configuration parameters.
        """
        self.fig = plt.figure(figsize=(10, 6))
        self.ax = self.fig.add_subplot(111)
        self.ax = self.plot_on_canvas(self.fig, self.ax, config)
        plt.tight_layout()
        plt.show()

    # ----- Helpers: evaluation & features -----
    def _evaluate_function(self, func_str, x):
        """Safely evaluates a mathematical expression string.

        A limited set of numpy functions and constants are made available
        to the expression.

        Args:
            func_str (str): The mathematical expression to evaluate, using 'x'
                as the independent variable.
            x (np.ndarray): The numpy array of x-values to evaluate against.

        Returns:
            np.ndarray: The result of the function evaluation (y-values).
        """
        func_str = str(func_str).replace('^', '**')
        safe = {
            'x': x,
            'sin': np.sin, 'cos': np.cos, 'tan': np.tan,
            'exp': np.exp, 'log': np.log, 'log10': np.log10,
            'sqrt': np.sqrt, 'abs': np.abs,
            'pi': np.pi, 'e': np.e
        }
        try:
            return eval(func_str, {"__builtins__": None}, safe)
        except Exception as e:
            print(f"Error evaluating '{func_str}': {e}")
            return np.zeros_like(x)

    def _central_diff(self, y, h):
        """Calculates the derivative using the central difference method.

        Args:
            y (np.ndarray): The array of y-values.
            h (float): The step size (delta x).

        Returns:
            np.ndarray: The calculated derivative, same length as y.
        """
        dy = np.empty_like(y)
        dy[1:-1] = (y[2:] - y[:-2]) / (2*h)
        dy[0] = (y[1] - y[0]) / h
        dy[-1] = (y[-1] - y[-2]) / h
        return dy

    def _zero_crossings(self, arr):
        """Finds the indices just before a zero-crossing in an array.

        A zero-crossing is detected where the sign of the array values changes.

        Args:
            arr (np.ndarray): The input array.

        Returns:
            np.ndarray: An array of indices.
        """
        s = np.sign(arr)
        return np.where(np.diff(s) != 0)[0]

    def _interp_zero(self, x0, x1, y0, y1):
        """Finds the x-value of a zero-crossing using linear interpolation.

        Args:
            x0 (float): The x-coordinate of the first point.
            x1 (float): The x-coordinate of the second point.
            y0 (float): The y-coordinate of the first point.
            y1 (float): The y-coordinate of the second point.

        Returns:
            float: The interpolated x-value where y is zero.
        """
        if y1 == y0:
            return (x0 + x1) / 2.0
        return x0 - y0 * (x1 - x0) / (y1 - y0)

    def _find_and_mark_features(self, x, y, h, config):
        """Finds and marks roots, extrema, and inflection points on the plot.

        Args:
            x (np.ndarray): The array of x-values.
            y (np.ndarray): The array of y-values.
            h (float): The step size (delta x).
            config (dict): The configuration dictionary, used to check which
                features to find (e.g., 'find_roots').
        """
        if config.get('find_roots'):
            idx = self._zero_crossings(y)
            xs = [self._interp_zero(x[i], x[i+1], y[i], y[i+1]) for i in idx]
            ys = np.zeros_like(xs, dtype=float)
            self.ax.scatter(xs, ys, marker='o', s=40, label='Roots')

        if config.get('find_extrema') or config.get('find_inflections'):
            dy = self._central_diff(y, h)

        if config.get('find_extrema'):
            idx = self._zero_crossings(dy)
            xs = [self._interp_zero(x[i], x[i+1], dy[i], dy[i+1]) for i in idx]
            ys = []
            for xi in xs:
                j = np.searchsorted(x, xi) - 1
                j = np.clip(j, 0, len(x)-2)
                ys.append(np.interp(xi, [x[j], x[j+1]], [y[j], y[j+1]]))
            self.ax.scatter(xs, ys, marker='^', s=50, label='Extrema')

        if config.get('find_inflections'):
            d2y = self._central_diff(self._central_diff(y, h), h)
            idx = self._zero_crossings(d2y)
            xs = [self._interp_zero(x[i], x[i+1], d2y[i], d2y[i+1]) for i in idx]
            ys = []
            for xi in xs:
                j = np.searchsorted(x, xi) - 1
                j = np.clip(j, 0, len(x)-2)
                ys.append(np.interp(xi, [x[j], x[j+1]], [y[j], y[j+1]]))
            self.ax.scatter(xs, ys, marker='s', s=45, label='Inflection')

    def _plot_data_file(self, config):
        """Overlays data from a text/CSV file onto the current plot.

        The data file is expected to have two columns (x, y). It tries to
        read with whitespace delimiters first, then falls back to commas.

        Args:
            config (dict): The configuration dictionary containing data file
                path and plotting options ('data_file', 'plot_type', etc.).
        """
        path = str(config.get('data_file'))
        if not os.path.isfile(path):
            print(f"Data file not found: {path}")
            return
        # Flexible delimiter: try whitespace first, then comma
        try:
            try:
                data = np.loadtxt(path, delimiter=None)
            except Exception:
                data = np.loadtxt(path, delimiter=',')
        except Exception as e:
            print(f"Could not read data file '{path}': {e}")
            return

        if data.ndim != 2 or data.shape[1] < 2:
            print("Data file must have at least two columns (x, y).")
            return

        x, y = data[:, 0], data[:, 1]
        plot_type = str(config.get('plot_type', 'scatter')).lower()
        color = config.get('data_color', 'tab:gray')
        label = config.get('data_label', 'data')

        if plot_type == 'line':
            self.ax.plot(x, y, linewidth=2, label=label, color=color)
        else:
            self.ax.scatter(x, y, label=label, s=16, color=color)

    def _configure_plot(self, config):
        """Applies final cosmetic settings to the plot.

        Sets title, labels, grid, and legend based on the config.

        Args:
            config (dict): The configuration dictionary.
        """
        self.ax.set_xlabel(config.get('xlabel', 'x'))
        self.ax.set_ylabel(config.get('ylabel', 'y'))
        self.ax.set_title(config.get('title', 'Function Plot'))
        if config.get('grid', True):
            self.ax.grid(True, alpha=0.3)
        if config.get('legend', True):
            self.ax.legend()

    # ----- Parametric -----
    def _plot_parametric(self, config):
        """Plots a parametric curve defined by x(t) and y(t).

        Args:
            config (dict): The configuration dictionary containing the
                parametric expressions 'px' and 'py' and t-range.
        """
        px = str(config.get('px'))
        py = str(config.get('py'))
        t0 = float(config.get('tstart', 0))
        t1 = float(config.get('tend', 2*np.pi))
        ts = float(config.get('tstep', 0.01))
        if ts <= 0:
            raise ValueError("tstep must be positive.")
        t = np.arange(t0, t1 + ts, ts, dtype=float)

        safe = {
            't': t,
            'sin': np.sin, 'cos': np.cos, 'tan': np.tan,
            'exp': np.exp, 'log': np.log, 'log10': np.log10,
            'sqrt': np.sqrt, 'abs': np.abs,
            'pi': np.pi, 'e': np.e
        }
        try:
            x = eval(px.replace('^', '**'), {"__builtins__": None}, safe)
            y = eval(py.replace('^', '**'), {"__builtins__": None}, safe)
        except Exception as e:
            print(f"Error evaluating parametric functions: {e}")
            x = np.zeros_like(t); y = np.zeros_like(t)

        self.ax.plot(x, y, label='parametric')
        self.ax.set_xlabel(config.get('xlabel', 'x(t)'))
        self.ax.set_ylabel(config.get('ylabel', 'y(t)'))
        self.ax.set_title(config.get('title', 'Parametric Plot'))

    # ----- Bode -----
    def _maybe_plot_bode(self, figure, config):
        """Plots a Bode diagram if Bode-specific keys are in the config.

        If 'bode_num' and 'bode_den' are present, this method clears the
        figure, creates two subplots for magnitude and phase, and plots
        the Bode diagram.

        Args:
            figure (matplotlib.figure.Figure): The figure to draw on.
            config (dict): The configuration dictionary.

        Returns:
            bool: True if a Bode plot was drawn, False otherwise.
        """
        if config.get('bode_num') is None or config.get('bode_den') is None:
            return False

        num = np.array(config.get('bode_num'), dtype=float)  # descending powers
        den = np.array(config.get('bode_den'), dtype=float)
        wmin = float(config.get('wmin', 0.1))
        wmax = float(config.get('wmax', 100.0))
        npts = int(config.get('npoints', 400))
        w = np.logspace(np.log10(wmin), np.log10(wmax), npts)

        jw = 1j * w
        Num = np.polyval(num, jw)
        Den = np.polyval(den, jw)
        H = Num / Den
        mag = 20*np.log10(np.abs(H) + 1e-30)  # dB
        phase = np.angle(H, deg=True)         # degrees

        figure.clear()
        ax1 = figure.add_subplot(2, 1, 1)
        ax2 = figure.add_subplot(2, 1, 2)
        ax1.semilogx(w, mag)
        ax2.semilogx(w, phase)
        ax1.set_title(config.get('title', 'Bode Diagram'))
        ax1.set_ylabel('Magnitude (dB)')
        ax1.grid(True, which='both', alpha=0.3)
        ax2.set_xlabel('Frequency (rad/s)')
        ax2.set_ylabel('Phase (deg)')
        ax2.grid(True, which='both', alpha=0.3)
        return True

    # ----- 3D / Contour -----
    def _maybe_plot_surface_or_contour(self, figure, config):
        """Plots a 3D surface and/or contour plot if 'fxy' is in config.

        If 'fxy' is present, this method clears the figure and creates either
        one or two subplots to show the surface and/or contour plot.

        Args:
            figure (matplotlib.figure.Figure): The figure to draw on.
            config (dict): The configuration dictionary.

        Returns:
            bool: True if a 3D plot was drawn, False otherwise.
        """
        fxy = config.get('fxy')
        if fxy is None:
            return False

        x0 = float(config.get('xstart', -2))
        x1 = float(config.get('xend', 2))
        xs = float(config.get('xstep', 0.05))
        y0 = float(config.get('ystart', -2))
        y1 = float(config.get('yend', 2))
        ys = float(config.get('ystep', 0.05))
        X = np.arange(x0, x1 + xs, xs, dtype=float)
        Y = np.arange(y0, y1 + ys, ys, dtype=float)
        XX, YY = np.meshgrid(X, Y)

        safe = {
            'x': XX, 'y': YY,
            'sin': np.sin, 'cos': np.cos, 'tan': np.tan,
            'exp': np.exp, 'log': np.log, 'log10': np.log10,
            'sqrt': np.sqrt, 'abs': np.abs,
            'pi': np.pi, 'e': np.e
        }
        try:
            ZZ = eval(str(fxy).replace('^', '**'), {"__builtins__": None}, safe)
        except Exception as e:
            print(f"Error evaluating f(x,y): {e}")
            ZZ = np.zeros_like(XX)

        figure.clear()
        surface_flag = bool(config.get('surface', True))
        contour_flag = bool(config.get('contour', True))
        levels = int(config.get('levels', 20))

        if surface_flag:
            ax3d = figure.add_subplot(1, 2 if contour_flag else 1, 1, projection='3d')
            ax3d.plot_surface(XX, YY, ZZ, linewidth=0, antialiased=True)
            ax3d.set_title('Surface')
            ax3d.set_xlabel('x'); ax3d.set_ylabel('y'); ax3d.set_zlabel(config.get('zlabel', 'f(x,y)'))

        if contour_flag:
            axc = figure.add_subplot(1, 2, 2) if surface_flag else figure.add_subplot(1, 1, 1)
            cs = axc.contourf(XX, YY, ZZ, levels=levels)
            axc.set_title('Contour')
            axc.set_xlabel('x'); axc.set_ylabel('y')
            figure.colorbar(cs, ax=axc, shrink=0.8)
        return True
