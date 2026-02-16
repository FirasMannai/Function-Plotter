
# widgets.py
# Custom Tk widgets and helpers.

import tkinter as tk
from tkinter import ttk

class ScrollableFrame(ttk.Frame):
    """A ttk.Frame that scrolls vertically (for tall control panels)."""
    def __init__(self, parent, *args, **kwargs):
        """Initializes the scrollable frame.

        This sets up a Canvas widget with a vertical scrollbar. An inner Frame
        is placed inside the canvas, and all content should be added to the
        'self.inner' frame. The scrolling mechanism is configured to move
        the canvas view over the inner frame.

        Args:
            parent: The parent tkinter widget.
            *args, **kwargs: Additional arguments passed to the parent ttk.Frame.
        """
        super().__init__(parent, *args, **kwargs)
        self.canvas = tk.Canvas(self, borderwidth=0, highlightthickness=0)
        self.vsb = ttk.Scrollbar(self, orient="vertical", command=self.canvas.yview)
        self.canvas.configure(yscrollcommand=self.vsb.set)
        self.vsb.pack(side="right", fill="y")
        self.canvas.pack(side="left", fill="both", expand=True)

        self.inner = ttk.Frame(self.canvas)
        self.inner_id = self.canvas.create_window((0, 0), window=self.inner, anchor="nw")

        self.inner.bind("<Configure>", lambda e: self.canvas.configure(scrollregion=self.canvas.bbox("all")))
        self.canvas.bind("<Configure>", lambda e: self.canvas.itemconfigure(self.inner_id, width=e.width))

        # Mouse wheel binding was removed to prevent critical UI stability issues.
        # Please use the vertical scrollbar to navigate the controls.
