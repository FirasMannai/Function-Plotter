
# app.py
# Entry point for GUI and headless modes.


import sys
import tkinter as tk
from gui import FunctionPlotterGUI
from config_parser import ConfigParser
from plotter import FunctionPlotter

def main():
    """Entry point for the Function Plotter application.

    The application can run in two modes:
    - GUI mode: If no command-line arguments are provided, the Tkinter
      graphical user interface will be launched.
      Example: `python app.py`

    - Headless mode: If a path to a configuration file is provided as a
      command-line argument, the application will parse the file and
      directly generate a plot in a Matplotlib window without launching the GUI.
      Example: `python app.py config.txt`
    """
    if len(sys.argv) > 1:
        # Headless: parse a config file and open a Matplotlib window
        cfg_file = sys.argv[1]
        parser = ConfigParser()
        cfg = parser.parse_config(cfg_file)
        if cfg:
            FunctionPlotter().plot_from_config(cfg)
        else:
            print("Error: Could not parse configuration file")
    else:
        # GUI mode
        root = tk.Tk()
        FunctionPlotterGUI(root)
        root.mainloop()

if __name__ == "__main__":
    main()
