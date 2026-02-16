
# config_parser.py
# Simple "key: value" configuration file reader.

import ast
import re

class ConfigParser:
    """Parses simple 'key: value' config files into a dictionary.

    This class reads text-based configuration files, handling comments (#),
    and converts values to appropriate Python types like bool, int, float, or list.
    It also handles special cases like collecting 'f1', 'f2', ... into a list.
    """
    def __init__(self):
        """Initializes the parser."""
        # A set of all known/supported configuration keys.
        # This is for documentation purposes and potential future validation.
        self.valid_keys = {
            'fi', 'functions',
            'xstart', 'xend', 'xstep',
            'ystart', 'yend', 'ystep',
            'grid', 'legend', 'title', 'xlabel', 'ylabel', 'zlabel',
            'color', 'colors', 'linewidth', 'linewidths', 'linestyles',
            'data_file', 'plot_type', 'data_color', 'data_label',
            'derivative1', 'derivative2',
            'find_roots', 'find_extrema', 'find_inflections',
            # Parametric
            'px', 'py', 'tstart', 'tend', 'tstep',
            # Bode
            'bode_num', 'bode_den', 'wmin', 'wmax', 'npoints',
            # 3D
            'fxy', 'surface', 'contour', 'levels'
        }

    def parse_config(self, filename):
        """Parses a configuration file and returns a dictionary of parameters.

        Args:
            filename (str): The path to the configuration file.

        Returns:
            dict or None: A dictionary containing the parsed configuration
            parameters, or None if an error occurs.
        """
        config = {}
        try:
            with open(filename, 'r', encoding='utf-8') as f:
                lines = f.readlines()

            for line in lines:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                if ':' in line:
                    key, value = line.split(':', 1)
                    key, value = key.strip(), value.strip()
                    config[key] = self._convert_value(value)

            # Collect f1, f2, ... into functions list (ascending order)
            func_keys = sorted(
                [k for k in config.keys() if re.fullmatch(r"f\d+", k)],
                key=lambda s: int(s[1:])
            )
            if func_keys:
                config['functions'] = [config[k] for k in func_keys]

            # Map fi -> functions when not already set
            if 'fi' in config and 'functions' not in config:
                config['functions'] = [config['fi']]

            return config
        except Exception as e:
            print(f"Error reading configuration file: {e}")
            return None

    def _convert_value(self, value):
        """Converts a string value to a more specific type if possible.

        Handles boolean strings ('true', 'false'), numeric values, and list-like
        strings (e.g., '[1, 2, 3]').

        Args:
            value (str): The string value from the config file.

        Returns:
            bool, int, float, list, or str: The converted value.
        """
        if isinstance(value, str):
            low = value.lower()
            if low in ('true', 'yes', 'on'):
                return True
            if low in ('false', 'no', 'off'):
                return False

        # numeric
        try:
            return float(value) if '.' in str(value) else int(value)
        except (ValueError, TypeError):
            pass

        # Python literal list (e.g., "[1, 2, 3]")
        if isinstance(value, str) and value.startswith('[') and value.endswith(']'):
            try:
                return ast.literal_eval(value)
            except Exception:
                # Fallback for list of non-literals, like expressions or unquoted strings
                content = value.strip()[1:-1].strip()
                if content:
                    return [item.strip() for item in content.split(',')]
                return []

        return value
