#!/usr/bin/env python3
"""GUI for photometry transformation and AAVSO report generation."""

import tkinter as tk
from tkinter import ttk, filedialog, scrolledtext, messagebox
from pathlib import Path
import threading
from datetime import datetime
import matplotlib
matplotlib.use('TkAgg')
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import sys
import os

# Import the analysis functions from analyze_photometry
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analyze_photometry import (
    read_table, extract_target_position, query_aavso_comparisons,
    match_comparison_stars, compute_transformation, apply_transformation_to_targets
)


class PhotometryGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Photometry Transformation Tool")
        self.root.geometry("1200x800")
        
        # Dark mode colors
        self.bg_dark = '#2b2b2b'
        self.bg_medium = '#3c3f41'
        self.bg_light = '#4a4d50'
        self.fg_main = '#ffffff'
        self.fg_dim = '#bbbbbb'
        self.accent_blue = '#1e88e5'  # Brighter, more saturated blue
        self.accent_green = '#43a047'  # Brighter, more saturated green
        
        # Configure root window
        self.root.configure(bg=self.bg_dark)
        
        # Set up dark theme
        style = ttk.Style()
        style.theme_use('default')
        
        # Configure ttk styles for dark mode
        style.configure('TFrame', background=self.bg_dark)
        style.configure('TLabel', background=self.bg_dark, foreground=self.fg_main, font=('TkDefaultFont', 10))
        style.configure('TLabelframe', background=self.bg_dark, foreground=self.fg_main)
        style.configure('TLabelframe.Label', background=self.bg_dark, foreground=self.fg_main, font=('TkDefaultFont', 10, 'bold'))
        style.configure('TNotebook', background=self.bg_dark, borderwidth=0)
        style.configure('TNotebook.Tab', background=self.bg_medium, foreground=self.fg_main, padding=[20, 10])
        style.map('TNotebook.Tab', background=[('selected', self.bg_light)], foreground=[('selected', self.fg_main)])
        style.configure('TEntry', fieldbackground=self.bg_light, foreground=self.fg_main, borderwidth=1)
        style.configure('TButton', background=self.bg_medium, foreground=self.fg_main)
        style.map('TButton', background=[('active', self.bg_light)])
        style.configure('TCombobox', fieldbackground=self.bg_light, foreground=self.fg_main, background=self.bg_medium)
        style.map('TCombobox', fieldbackground=[('readonly', self.bg_light)], foreground=[('readonly', self.fg_main)])
        
        # Variables
        self.tbl_file = tk.StringVar()
        # Set default transform coefficients file
        default_coeffs = Path(__file__).parent / "transform_coeffs_short.txt"
        self.coeffs_file = tk.StringVar(value=str(default_coeffs) if default_coeffs.exists() else "")
        self.target_auid = tk.StringVar(value="000-BQH-515")
        self.output_file = tk.StringVar()
        self.obscode = tk.StringVar(value="KSCA")
        self.color_basis = tk.StringVar(value="vi")
        self.match_tolerance = tk.DoubleVar(value=2.0/60.0)  # 2 arcsec in arcmin
        self.aavso_radius = tk.DoubleVar(value=30.0)
        self.aavso_mag_limit = tk.DoubleVar(value=15.0)
        
        # Data storage
        self.rows = None
        self.transformation = None
        self.target_mags = None
        self.matches = None
        self.aavso_payload = None  # Store AAVSO response for chart ID
        
        self.create_widgets()
        
    def create_widgets(self):
        # Main container with notebook tabs
        notebook = ttk.Notebook(self.root)
        notebook.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # Tab 1: Configuration
        config_frame = ttk.Frame(notebook)
        notebook.add(config_frame, text="Configuration")
        self.create_config_tab(config_frame)
        
        # Tab 2: Time Series Plot
        plot_frame = ttk.Frame(notebook)
        notebook.add(plot_frame, text="Time Series Plot")
        self.create_plot_tab(plot_frame)
        
        # Tab 3: Results
        results_frame = ttk.Frame(notebook)
        notebook.add(results_frame, text="Results")
        self.create_results_tab(results_frame)
        
    def create_config_tab(self, parent):
        # File selection section
        file_frame = ttk.LabelFrame(parent, text="Input Files", padding=10)
        file_frame.pack(fill=tk.X, padx=10, pady=5)
        
        # Photometry table
        row = 0
        ttk.Label(file_frame, text="Photometry Table (.tbl):", font=('TkDefaultFont', 10)).grid(row=row, column=0, sticky=tk.W, pady=5, padx=5)
        ttk.Entry(file_frame, textvariable=self.tbl_file, width=50).grid(row=row, column=1, padx=5)
        ttk.Button(file_frame, text="Browse", command=self.browse_tbl, width=10).grid(row=row, column=2, padx=5)
        
        # Transform coefficients (mandatory)
        row += 1
        ttk.Label(file_frame, text="Transform Coefficients:", font=('TkDefaultFont', 10)).grid(row=row, column=0, sticky=tk.W, pady=5, padx=5)
        ttk.Entry(file_frame, textvariable=self.coeffs_file, width=50).grid(row=row, column=1, padx=5)
        ttk.Button(file_frame, text="Browse", command=self.browse_coeffs, width=10).grid(row=row, column=2, padx=5)
        
        # Target configuration
        target_frame = ttk.LabelFrame(parent, text="Target Configuration", padding=10)
        target_frame.pack(fill=tk.X, padx=10, pady=5)
        
        ttk.Label(target_frame, text="Target AUID:", font=('TkDefaultFont', 10)).grid(row=0, column=0, sticky=tk.W, pady=5, padx=5)
        ttk.Entry(target_frame, textvariable=self.target_auid, width=20, font=('TkDefaultFont', 10)).grid(row=0, column=1, sticky=tk.W, padx=5)
        
        ttk.Label(target_frame, text="Observer Code:", font=('TkDefaultFont', 10)).grid(row=1, column=0, sticky=tk.W, pady=5, padx=5)
        ttk.Entry(target_frame, textvariable=self.obscode, width=20, font=('TkDefaultFont', 10)).grid(row=1, column=1, sticky=tk.W, padx=5)
        
        ttk.Label(target_frame, text="Color Basis:", font=('TkDefaultFont', 10)).grid(row=2, column=0, sticky=tk.W, pady=5, padx=5)
        color_combo = ttk.Combobox(target_frame, textvariable=self.color_basis, 
                                    values=["bv", "vr", "vi"], state="readonly", width=18, font=('TkDefaultFont', 10))
        color_combo.grid(row=2, column=1, sticky=tk.W, padx=5)
        
        # AAVSO parameters
        aavso_frame = ttk.LabelFrame(parent, text="AAVSO Query Parameters", padding=10)
        aavso_frame.pack(fill=tk.X, padx=10, pady=5)
        
        ttk.Label(aavso_frame, text="Search Radius (arcmin):", font=('TkDefaultFont', 10)).grid(row=0, column=0, sticky=tk.W, pady=5, padx=5)
        ttk.Entry(aavso_frame, textvariable=self.aavso_radius, width=20, font=('TkDefaultFont', 10)).grid(row=0, column=1, sticky=tk.W, padx=5)
        
        ttk.Label(aavso_frame, text="Magnitude Limit:", font=('TkDefaultFont', 10)).grid(row=1, column=0, sticky=tk.W, pady=5, padx=5)
        ttk.Entry(aavso_frame, textvariable=self.aavso_mag_limit, width=20, font=('TkDefaultFont', 10)).grid(row=1, column=1, sticky=tk.W, padx=5)
        
        ttk.Label(aavso_frame, text="Match Tolerance (arcsec):", font=('TkDefaultFont', 10)).grid(row=2, column=0, sticky=tk.W, pady=5, padx=5)
        tolerance_entry = ttk.Entry(aavso_frame, width=20, font=('TkDefaultFont', 10))
        tolerance_entry.insert(0, "2.0")
        tolerance_entry.grid(row=2, column=1, sticky=tk.W, padx=5)
        tolerance_entry.bind('<FocusOut>', lambda e: self.match_tolerance.set(float(tolerance_entry.get())/60.0))
        
        # Output section
        output_frame = ttk.LabelFrame(parent, text="Output", padding=10)
        output_frame.pack(fill=tk.X, padx=10, pady=5)
        
        ttk.Label(output_frame, text="AAVSO Report Output:", font=('TkDefaultFont', 10)).grid(row=0, column=0, sticky=tk.W, pady=5, padx=5)
        ttk.Entry(output_frame, textvariable=self.output_file, width=50).grid(row=0, column=1, padx=5)
        ttk.Button(output_frame, text="Browse", command=self.browse_output, width=10).grid(row=0, column=2, padx=5)
        
        # Action buttons
        button_frame = ttk.Frame(parent)
        button_frame.pack(fill=tk.X, padx=10, pady=10)
        
        self.process_btn = tk.Button(button_frame, text="Process Transformation", 
                                      command=self.process_transformation,
                                      bg=self.accent_blue, fg='white', font=('TkDefaultFont', 11, 'bold'),
                                      padx=25, pady=12, relief=tk.RAISED, borderwidth=2,
                                      activebackground='#1565c0', activeforeground='white',
                                      cursor='hand2', highlightthickness=1, highlightbackground=self.accent_blue)
        self.process_btn.pack(side=tk.LEFT, padx=5)
        
        self.export_btn = tk.Button(button_frame, text="Export AAVSO Report", 
                                     command=self.export_aavso_report,
                                     bg=self.accent_green, fg='white', font=('TkDefaultFont', 11, 'bold'),
                                     padx=25, pady=12, relief=tk.RAISED, borderwidth=2,
                                     activebackground='#2e7d32', activeforeground='white',
                                     cursor='hand2', highlightthickness=1, highlightbackground=self.accent_green,
                                     state=tk.DISABLED)
        self.export_btn.pack(side=tk.LEFT, padx=5)
        
        # Status area
        status_frame = ttk.LabelFrame(parent, text="Status", padding=10)
        status_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)
        
        self.status_text = scrolledtext.ScrolledText(status_frame, height=10, wrap=tk.WORD,
                                                       bg=self.bg_light, fg=self.fg_main,
                                                       insertbackground=self.fg_main,
                                                       font=('Courier', 10))
        self.status_text.pack(fill=tk.BOTH, expand=True)
        
        # Initial status message
        self.log("Photometry Transformation Tool - Ready")
        self.log("1. Select a photometry table (.tbl file)")
        self.log("2. Configure target AUID and observer code")
        self.log("3. Verify transformation coefficients file is correct")
        self.log("4. Click 'Process Transformation' to begin")
        self.log("-" * 60)
        
    def create_plot_tab(self, parent):
        # Create matplotlib figure with dark theme
        self.figure = Figure(figsize=(10, 6), dpi=100, facecolor=self.bg_dark)
        self.ax = self.figure.add_subplot(111, facecolor=self.bg_medium)
        self.ax.set_xlabel('Observation Number', color=self.fg_main)
        self.ax.set_ylabel('V Magnitude', color=self.fg_main)
        self.ax.set_title('Time Series: T1 (Target) and T2 (Check Star)', color=self.fg_main)
        self.ax.grid(True, alpha=0.2, color=self.fg_dim)
        self.ax.tick_params(colors=self.fg_main)
        self.ax.spines['bottom'].set_color(self.fg_dim)
        self.ax.spines['top'].set_color(self.fg_dim)
        self.ax.spines['left'].set_color(self.fg_dim)
        self.ax.spines['right'].set_color(self.fg_dim)
        self.ax.invert_yaxis()  # Magnitudes are inverted
        
        # Create canvas
        self.canvas = FigureCanvasTkAgg(self.figure, parent)
        self.canvas.get_tk_widget().configure(bg=self.bg_dark)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        # Toolbar
        toolbar_frame = ttk.Frame(parent)
        toolbar_frame.pack(fill=tk.X, padx=10, pady=5)
        
        tk.Button(toolbar_frame, text="Refresh Plot", command=self.update_plot,
                  bg=self.bg_medium, fg=self.fg_main, font=('TkDefaultFont', 9),
                  padx=10, pady=5, relief=tk.FLAT, borderwidth=0,
                  activebackground=self.bg_light, activeforeground=self.fg_main,
                  cursor='hand2').pack(side=tk.LEFT, padx=5)
        
    def create_results_tab(self, parent):
        results_text_frame = ttk.Frame(parent)
        results_text_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        self.results_text = scrolledtext.ScrolledText(results_text_frame, wrap=tk.WORD, 
                                                       font=("Courier", 10),
                                                       bg=self.bg_light, fg=self.fg_main,
                                                       insertbackground=self.fg_main)
        self.results_text.pack(fill=tk.BOTH, expand=True)
        
    def browse_tbl(self):
        filename = filedialog.askopenfilename(
            parent=self.root,
            title="Select Photometry Table",
            filetypes=[("Table files", "*.tbl"), ("All files", "*.*")]
        )
        if filename:
            self.tbl_file.set(filename)
            # Auto-generate output filename
            base = Path(filename).stem
            out_path = Path(filename).parent / f"AAVSO_AID_Report_{self.obscode.get()}_{self.target_auid.get().replace('-', '')}_{datetime.now().strftime('%Y%m%d-%Hh%Mm%Ss')}.txt"
            self.output_file.set(str(out_path))
            self.log(f"Loaded: {filename}")
            
    def browse_coeffs(self):
        filename = filedialog.askopenfilename(
            parent=self.root,
            title="Select Transformation Coefficients",
            filetypes=[("Text files", "*.txt"), ("INI files", "*.ini"), ("All files", "*.*")]
        )
        if filename:
            self.coeffs_file.set(filename)
            self.log(f"Loaded coefficients: {filename}")
            
    def browse_output(self):
        filename = filedialog.asksaveasfilename(
            parent=self.root,
            title="Save AAVSO Report As",
            defaultextension=".txt",
            filetypes=[("Text files", "*.txt"), ("All files", "*.*")]
        )
        if filename:
            self.output_file.set(filename)
            
    def log(self, message):
        """Add message to status log."""
        self.status_text.insert(tk.END, message + "\n")
        self.status_text.see(tk.END)
        self.root.update_idletasks()
        
    def process_transformation(self):
        """Main processing function."""
        if not self.tbl_file.get():
            messagebox.showerror("Error", "Please select a photometry table file.")
            return
        
        if not self.coeffs_file.get():
            messagebox.showerror("Error", "Please select a transformation coefficients file.")
            return
        
        # Check if coefficients file exists
        coeffs_path = Path(self.coeffs_file.get())
        if not coeffs_path.exists():
            messagebox.showerror("Error", f"Transformation coefficients file not found:\n{coeffs_path}")
            return
            
        # Disable button during processing
        self.process_btn.config(state=tk.DISABLED)
        self.status_text.delete(1.0, tk.END)
        
        # Run in thread to avoid freezing GUI
        thread = threading.Thread(target=self._do_processing, daemon=True)
        thread.start()
        
    def _do_processing(self):
        """Worker thread for processing."""
        try:
            self.log("=" * 60)
            self.log("Starting photometry transformation...")
            self.log("=" * 60)
            
            # Read table
            self.log(f"\nReading table: {self.tbl_file.get()}")
            table_path = Path(self.tbl_file.get())
            header, self.rows = read_table(table_path)
            self.log(f"  Found {len(self.rows)} observations")
            
            # Extract target position
            self.log("\nExtracting target position...")
            ra_deg, dec_deg = extract_target_position(self.rows)
            self.log(f"  RA: {ra_deg:.6f}°")
            self.log(f"  Dec: {dec_deg:.6f}°")
            
            # Query AAVSO
            self.log("\nQuerying AAVSO for comparison stars...")
            payload = query_aavso_comparisons(
                ra_deg, dec_deg,
                radius_arcmin=self.aavso_radius.get(),
                mag_limit=self.aavso_mag_limit.get()
            )
            
            # Store payload for later use
            self.aavso_payload = payload
            
            aavso_stars = payload.get("photometry", [])
            self.log(f"  Retrieved {len(aavso_stars)} comparison stars")
            
            # Match comparison stars
            self.log("\nMatching comparison stars...")
            self.matches = match_comparison_stars(
                self.rows, aavso_stars, ra_deg, dec_deg,
                match_tolerance_arcmin=self.match_tolerance.get()
            )
            
            if not self.matches:
                self.log("  ERROR: No comparison stars matched!")
                self.root.after(0, lambda: self.process_btn.config(state=tk.NORMAL))
                return
                
            self.log(f"  Matched {len(self.matches)} stars:")
            for label, match_data in sorted(self.matches.items()):
                aavso_label = match_data.get("label", "?")
                sep = match_data.get("separation_arcmin", 0) * 60  # Convert to arcsec
                self.log(f"    {label} -> AAVSO {aavso_label} (sep={sep:.2f}\")")
            
            # Compute transformation
            self.log(f"\nComputing transformation (color basis: {self.color_basis.get().upper()})...")
            self.transformation = compute_transformation(
                self.rows, self.matches, 
                color_basis=self.color_basis.get()
            )
            
            zp = self.transformation["zero_point"]
            k = self.transformation["color_coefficient"]
            rms = self.transformation["rms"]
            n_stars = self.transformation["num_stars"]
            n_obs = self.transformation["num_obs"]
            
            self.log(f"  Zero point (ZP):        {zp:+.4f} mag")
            self.log(f"  Color coefficient (k):  {k:+.4f} mag/{self.color_basis.get().upper()}")
            self.log(f"  RMS residual:           {rms:.4f} mag")
            self.log(f"  Number of comp stars:   {n_stars}")
            self.log(f"  Number of observations: {n_obs}")
            
            # Apply transformation to targets
            self.log("\nApplying transformation to T1 and T2...")
            self.target_mags = apply_transformation_to_targets(
                self.rows, self.transformation, aavso_stars, self.matches,
                color_basis=self.color_basis.get()
            )
            
            # Display results
            for target in ["T1", "T2"]:
                if self.target_mags[target]:
                    mags = [mag for _, mag in self.target_mags[target]]
                    import statistics
                    mean_mag = statistics.fmean(mags)
                    std_dev = statistics.stdev(mags) if len(mags) > 1 else 0.0
                    self.log(f"\n{target} Results:")
                    self.log(f"  Mean V:  {mean_mag:.4f} ± {std_dev:.4f} mag")
                    self.log(f"  N obs:   {len(mags)}")
            
            self.log("\n" + "=" * 60)
            self.log("Processing complete!")
            self.log("=" * 60)
            
            # Update UI
            self.root.after(0, self.update_plot)
            self.root.after(0, self.display_results)
            self.root.after(0, lambda: self.export_btn.config(state=tk.NORMAL))
            
        except Exception as e:
            self.log(f"\nERROR: {str(e)}")
            import traceback
            self.log(traceback.format_exc())
            
        finally:
            self.root.after(0, lambda: self.process_btn.config(state=tk.NORMAL))
            
    def update_plot(self):
        """Update the time series plot."""
        if self.target_mags is None:
            return
            
        self.ax.clear()
        
        # Reapply dark theme after clear
        self.ax.set_facecolor(self.bg_medium)
        self.ax.tick_params(colors=self.fg_main)
        self.ax.spines['bottom'].set_color(self.fg_dim)
        self.ax.spines['top'].set_color(self.fg_dim)
        self.ax.spines['left'].set_color(self.fg_dim)
        self.ax.spines['right'].set_color(self.fg_dim)
        
        # Plot T1 (blue)
        if self.target_mags["T1"]:
            obs_nums_t1 = [idx + 1 for idx, _ in self.target_mags["T1"]]
            mags_t1 = [mag for _, mag in self.target_mags["T1"]]
            self.ax.plot(obs_nums_t1, mags_t1, color='#4fc3f7', linewidth=1.2, label='T1 (Target)', alpha=0.9)
            self.ax.scatter(obs_nums_t1, mags_t1, c='#4fc3f7', s=15, alpha=0.6, edgecolors='#29b6f6', linewidth=0.5)
        
        # Plot T2 (red/orange)
        if self.target_mags["T2"]:
            obs_nums_t2 = [idx + 1 for idx, _ in self.target_mags["T2"]]
            mags_t2 = [mag for _, mag in self.target_mags["T2"]]
            self.ax.plot(obs_nums_t2, mags_t2, color='#ff7043', linewidth=1.2, label='T2 (Check)', alpha=0.9)
            self.ax.scatter(obs_nums_t2, mags_t2, c='#ff7043', s=15, alpha=0.6, edgecolors='#ff5722', linewidth=0.5)
        
        self.ax.set_xlabel('Observation Number', color=self.fg_main, fontsize=11)
        self.ax.set_ylabel('V Magnitude', color=self.fg_main, fontsize=11)
        self.ax.set_title('Time Series: T1 (Target) and T2 (Check Star)', color=self.fg_main, fontsize=12, pad=15)
        legend = self.ax.legend(facecolor=self.bg_light, edgecolor=self.fg_dim, framealpha=0.9)
        for text in legend.get_texts():
            text.set_color(self.fg_main)
        self.ax.grid(True, alpha=0.2, color=self.fg_dim)
        self.ax.invert_yaxis()
        
        self.canvas.draw()
        
    def display_results(self):
        """Display detailed results in the results tab."""
        if self.target_mags is None:
            return
            
        self.results_text.delete(1.0, tk.END)
        
        # Header
        self.results_text.insert(tk.END, "=" * 80 + "\n")
        self.results_text.insert(tk.END, "PHOTOMETRY TRANSFORMATION RESULTS\n")
        self.results_text.insert(tk.END, "=" * 80 + "\n\n")
        
        # Transformation parameters
        if self.transformation:
            self.results_text.insert(tk.END, "Transformation Parameters:\n")
            self.results_text.insert(tk.END, f"  Zero Point:          {self.transformation['zero_point']:+.4f} mag\n")
            self.results_text.insert(tk.END, f"  Color Coefficient:   {self.transformation['color_coefficient']:+.4f} mag/{self.color_basis.get().upper()}\n")
            self.results_text.insert(tk.END, f"  RMS Residual:        {self.transformation['rms']:.4f} mag\n")
            self.results_text.insert(tk.END, f"  Comparison Stars:    {self.transformation['num_stars']}\n")
            self.results_text.insert(tk.END, f"  Observations Used:   {self.transformation['num_obs']}\n\n")
        
        # Individual measurements
        for target in ["T1", "T2"]:
            if self.target_mags[target]:
                self.results_text.insert(tk.END, f"\n{target} Measurements (first 50):\n")
                self.results_text.insert(tk.END, f"  {'Obs#':>5s}  {'V_std':>8s}\n")
                self.results_text.insert(tk.END, "  " + "-" * 16 + "\n")
                
                for i, (row_idx, v_std) in enumerate(self.target_mags[target][:50]):
                    self.results_text.insert(tk.END, f"  {row_idx+1:5d}  {v_std:8.4f}\n")
                    
                if len(self.target_mags[target]) > 50:
                    self.results_text.insert(tk.END, f"  ... ({len(self.target_mags[target]) - 50} more observations)\n")
                    
    def export_aavso_report(self):
        """Export results in AAVSO extended format."""
        if not self.target_mags or not self.output_file.get():
            messagebox.showerror("Error", "No results to export or output file not specified.")
            return
        
        if not self.matches:
            messagebox.showerror("Error", "No star matching data available. Run transformation first.")
            return
            
        try:
            output_path = Path(self.output_file.get())
            
            self.log(f"\nPreparing AAVSO report export...")
            self.log(f"  T1 observations: {len(self.target_mags.get('T1', []))}")
            self.log(f"  T2 observations: {len(self.target_mags.get('T2', []))}")
            self.log(f"  Matched stars: {list(self.matches.keys())}")
            
            # Determine CNAME and CMAG based on number of comparison stars
            num_comp_stars = len([k for k in self.matches.keys() if k.startswith('C')])
            
            if num_comp_stars > 1:
                cname = "ENSEMBLE"
                cmag = "na"
            elif num_comp_stars == 1:
                # Get the single comparison star
                comp_label = [k for k in self.matches.keys() if k.startswith('C')][0]
                comp_match = self.matches[comp_label]
                cname = comp_match.get("label", "na")
                # Get V magnitude from bands
                cmag = "na"
                for band in comp_match.get("bands", []):
                    if isinstance(band, dict) and band.get("band") == "V":
                        cmag = f"{band.get('mag', 0.0):.3f}"
                        break
            else:
                cname = "na"
                cmag = "na"
            
            self.log(f"  CNAME: {cname}, CMAG: {cmag}")
            
            # Get T2 (check star) name from matches
            kname = "na"
            if "T2" in self.matches:
                kname = self.matches["T2"].get("label", "na")
                self.log(f"  KNAME (check star): {kname}")
            else:
                self.log(f"  WARNING: T2 not in matches")
            
            # Get T2 (check star) transformed magnitude
            kmag = "na"
            if self.target_mags.get("T2"):
                import statistics
                t2_mags = [mag for _, mag in self.target_mags["T2"]]
                kmag = f"{statistics.fmean(t2_mags):.3f}"
            
            # Get chart ID from AAVSO payload
            chart_id = "na"
            if self.aavso_payload:
                chart_id = self.aavso_payload.get("chartid", "na")
            
            self.log(f"  Chart ID: {chart_id}")
            
            # Determine which date column is available and will be used
            date_column_used = None  # Track which column we're actually using
            if self.rows and len(self.rows) > 0:
                sample_row = self.rows[0]
                for col in ["HJD_UTC", "HJD", "JD", "hjd_utc", "hjd", "jd"]:
                    if col in sample_row:
                        try:
                            float(sample_row[col])
                            date_column_used = col
                            break
                        except (ValueError, TypeError):
                            continue
            
            # For #DATE= header: use "HJD" if column is HJD or HJD_UTC, otherwise use actual column name
            if date_column_used and date_column_used.upper() in ["HJD", "HJD_UTC"]:
                date_header = "HJD"
            elif date_column_used:
                date_header = date_column_used
            else:
                date_header = "HJD"  # default fallback
            
            self.log(f"  Using date column: {date_column_used} (header: {date_header})")
            
            rows_written = 0
            with output_path.open("w") as f:
                # Write header
                f.write("#TYPE=Extended\n")
                f.write(f"#OBSCODE={self.obscode.get()}\n")
                f.write("#SOFTWARE=MyPhot Photometry Transformation Tool v1.0\n")
                f.write("#DELIM=,\n")
                f.write(f"#DATE={date_header}\n")
                f.write("#OBSTYPE=CCD\n")
                f.write("#NAME,DATE,MAG,MERR,FILT,TRANS,MTYPE,CNAME,CMAG,KNAME,KMAG,AMASS,GROUP,CHART,NOTES\n")
                
                # Write T1 observations
                if not self.target_mags.get("T1"):
                    self.log("  ERROR: No T1 data available!")
                    messagebox.showerror("Error", "No T1 target magnitudes available for export.")
                    return
                
                for row_idx, v_std in self.target_mags["T1"]:
                    if row_idx >= len(self.rows):
                        self.log(f"  WARNING: Row index {row_idx} out of range")
                        continue
                        
                    row = self.rows[row_idx]
                    
                    # Extract HJD_UTC (preferred) or fall back to other HJD columns
                    hjd = None
                    for hjd_col in ["HJD_UTC", "HJD", "JD", "hjd_utc", "hjd", "jd"]:
                        if hjd_col in row:
                            try:
                                hjd = float(row[hjd_col])
                                break
                            except (ValueError, TypeError):
                                continue
                    
                    if hjd is None:
                        self.log(f"  WARNING: No HJD found for row {row_idx}, skipping")
                        continue
                    
                    # Get individual T2 magnitude for this observation
                    kmag_obs = kmag  # Use mean T2 mag as default
                    if self.target_mags.get("T2"):
                        for t2_idx, t2_mag in self.target_mags["T2"]:
                            if t2_idx == row_idx:
                                kmag_obs = f"{t2_mag:.3f}"
                                break
                    
                    # Magnitude error: combine photometric error with transformation RMS in quadrature
                    # Get photometric error from Source_AMag_Err_T1
                    phot_err = 0.0
                    err_col_found = None
                    for err_col in ["Source_AMag_Err_T1", "Source_Error_T1", "source_amag_err_t1"]:
                        if err_col in row:
                            try:
                                phot_err = float(row[err_col])
                                err_col_found = err_col
                                break
                            except (ValueError, TypeError):
                                continue
                    
                    # Transformation RMS error
                    trans_rms = self.transformation['rms'] if self.transformation else 0.010
                    
                    # Combine in quadrature: total_err = sqrt(phot_err^2 + trans_rms^2)
                    import math
                    total_err = math.sqrt(phot_err**2 + trans_rms**2)
                    merr = f"{total_err:.3f}"
                    
                    # Debug logging for first few rows
                    if rows_written < 3:
                        self.log(f"  Row {row_idx}: phot_err={phot_err:.6f} (from {err_col_found}), trans_rms={trans_rms:.4f}, total={total_err:.6f} → {merr}")
                    
                    # Get airmass
                    amass = 1.0
                    for amass_col in ["AIRMASS", "Airmass", "airmass", "AirMass"]:
                        if amass_col in row:
                            try:
                                amass = float(row[amass_col])
                                break
                            except (ValueError, TypeError):
                                continue
                    
                    # Write data line
                    f.write(f"{self.target_auid.get()},")
                    f.write(f"{hjd:.6f},")
                    f.write(f"{v_std:.3f},")
                    f.write(f"{merr},")
                    f.write(f"V,")
                    f.write(f"YES,")
                    f.write(f"STD,")
                    f.write(f"{cname},")
                    f.write(f"{cmag},")
                    f.write(f"{kname},")
                    f.write(f"{kmag_obs},")
                    f.write(f"{amass:.6f},")
                    f.write(f"1,")
                    f.write(f"{chart_id},")
                    f.write(f"\n")  # Notes blank as specified
                    
                    rows_written += 1
            
            self.log(f"\nAAVSO report exported to: {output_path}")
            self.log(f"  Total observations written: {rows_written}")
            self.log(f"  Comparison stars: {cname}")
            self.log(f"  Check star: {kname} = {kmag}")
            self.log(f"  Chart ID: {chart_id}")
            
            if rows_written == 0:
                messagebox.showwarning("Warning", f"AAVSO report created but no data rows were written!\nCheck the status log for details.\n\n{output_path}")
            else:
                messagebox.showinfo("Success", f"AAVSO report exported successfully!\n\n{rows_written} observations written to:\n{output_path}")
            
        except Exception as e:
            messagebox.showerror("Error", f"Failed to export AAVSO report:\n{str(e)}")
            self.log(f"\nERROR exporting report: {str(e)}")
            import traceback
            self.log(traceback.format_exc())


def main():
    root = tk.Tk()
    app = PhotometryGUI(root)
    root.mainloop()


if __name__ == "__main__":
    main()
