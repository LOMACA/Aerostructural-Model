import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext
from airfoil_aerodynamic_data import AirfoilData

class AirfoilDataDisplay(ttk.Frame):
    def __init__(self, container, **kwargs):
        super().__init__(container, **kwargs)
        self.airfoil_data = AirfoilData()

        title_label = ttk.Label(self, text="Airfoil Data Import", font=("Arial", 16, "bold"))
        title_label.grid(row=0, column=0, columnspan=4, padx=10, pady=10)

        self.angle_start_label = ttk.Label(self, text="Enter Start Angle of Attack:")
        self.angle_start_label.grid(row=1, column=0, padx=5, pady=5)
        self.angle_start_entry = ttk.Entry(self)
        self.angle_start_entry.grid(row=1, column=1, padx=5, pady=5)

        self.angle_end_label = ttk.Label(self, text="Enter End Angle of Attack:")
        self.angle_end_label.grid(row=2, column=0, padx=5, pady=5)
        self.angle_end_entry = ttk.Entry(self)
        self.angle_end_entry.grid(row=2, column=1, padx=5, pady=5)

        ttk.Label(self, text="Tip Airfoil File:").grid(row=3, column=0, padx=5, pady=5)
        self.tip_file_entry = ttk.Entry(self, width=40)
        self.tip_file_entry.grid(row=3, column=1, padx=5, pady=5)
        ttk.Button(self, text="Browse", command=lambda: self.browse_file(self.tip_file_entry)).grid(row=3, column=2, padx=5)

        ttk.Label(self, text="Root Airfoil File:").grid(row=4, column=0, padx=5, pady=5)
        self.root_file_entry = ttk.Entry(self, width=40)
        self.root_file_entry.grid(row=4, column=1, padx=5, pady=5)
        ttk.Button(self, text="Browse", command=lambda: self.browse_file(self.root_file_entry)).grid(row=4, column=2, padx=5)

        title_label = ttk.Label(self, text="If no data is available: ", font=("Arial", 16, "bold"))
        title_label.grid(row=0, column=4, columnspan=4, padx=10, pady=10)

        ttk.Label(self, text="Tip Reynolds Number:").grid(row=1, column=4, padx=5, pady=5)
        self.tip_reynolds_entry = ttk.Entry(self)
        self.tip_reynolds_entry.grid(row=1, column=5, padx=5, pady=5)

        ttk.Label(self, text="Root Reynolds Number:").grid(row=2, column=4, padx=5, pady=5)
        self.root_reynolds_entry = ttk.Entry(self)
        self.root_reynolds_entry.grid(row=2, column=5, padx=5, pady=5)

        ttk.Label(self, text="Tip Airfoil File:").grid(row=3, column=4, padx=5, pady=5)
        self.tip_file_XFoil_entry = ttk.Entry(self, width=40)
        self.tip_file_XFoil_entry.grid(row=3, column=5, padx=5, pady=5)
        ttk.Button(self, text="Browse", command=lambda: self.browse_file(self.tip_file_XFoil_entry)).grid(row=3, column=6, padx=5)

        ttk.Label(self, text="Root Airfoil File:").grid(row=4, column=4, padx=5, pady=5)
        self.root_file_XFoil_entry = ttk.Entry(self, width=40)
        self.root_file_XFoil_entry.grid(row=4, column=5, padx=5, pady=5)
        ttk.Button(self, text="Browse", command=lambda: self.browse_file(self.root_file_XFoil_entry)).grid(row=4, column=6, padx=5)

        ttk.Button(self, text="Import Data", command=self.import_data).grid(row=5, column=0, columnspan=2, padx=5, pady=5)
        ttk.Button(self, text="Run XFoil", command=self.run_xfoil).grid(row=5, column=4, columnspan=2, padx=5, pady=5)

        self.status_label = ttk.Label(self, text="Status: Waiting for files...")
        self.status_label.grid(row=5, column=1, columnspan=4, padx=5, pady=5)

        self.data_display = scrolledtext.ScrolledText(self, width=40, height=10)
        self.data_display.grid(row=6, column=0, columnspan=4, padx=5, pady=5)

        self.status_label_Xfoil = ttk.Label(self, text="Waiting for XFoil...")
        self.status_label_Xfoil.grid(row=5, column=5, columnspan=4, padx=5, pady=5)

        self.data_display_XFoil = scrolledtext.ScrolledText(self, width=40, height=10)
        self.data_display_XFoil.grid(row=6, column=4, columnspan=4, padx=5, pady=5)

    def browse_file(self, entry_field):
        file_path = filedialog.askopenfilename()
        if file_path:
            entry_field.delete(0, tk.END)
            entry_field.insert(0, file_path)

    def import_data(self):
        tip_file_path = self.tip_file_entry.get()
        root_file_path = self.root_file_entry.get()

        try:
            angle_start = float(self.angle_start_entry.get())
            angle_end = float(self.angle_end_entry.get())
        except ValueError:
            messagebox.showerror("Error", "Invalid angle values entered. Please enter numeric values.")
            return

        tip_data = self.airfoil_data.data_extraction(tip_file_path, angle_start, angle_end)
        root_data = self.airfoil_data.data_extraction(root_file_path, angle_start, angle_end)

        if not tip_data:
            messagebox.showerror("Error",
                                 "No data found for the specified range of angles of attack for the tip airfoil.")
            return
        if not root_data:
            messagebox.showerror("Error",
                                 "No data found for the specified range of angles of attack for the root airfoil.")
            return

        if 'CL' not in tip_data[0] or 'CD' not in tip_data[0]:
            messagebox.showerror("Error", "Data format error in the tip airfoil.")
            return
        if 'CL' not in root_data[0] or 'CD' not in root_data[0]:
            messagebox.showerror("Error", "Data format error in the root airfoil.")
            return


        self.tip_cl0 = tip_data[0]['CL']
        self.tip_cd0 = tip_data[0]['CD']
        self.root_cl0 = root_data[0]['CL']
        self.root_cd0 = root_data[0]['CD']

        self.display_data(tip_data, "Tip Airfoil")
        self.display_data(root_data, "Root Airfoil")

    def run_xfoil(self):
        # This function would run XFoil with the parameters provided
        pass

    def display_data(self, data, airfoil_description):
        self.data_display.insert(tk.END, f"{airfoil_description} Data:\n")
        if data:
            for entry in data:
                self.data_display.insert(tk.END,
                                         f"Alpha: {entry['alpha']}, cl: {entry['CL']}, cd: {entry['CD']}, cdp: {entry['CDp']}, cm: {entry['CM']}, Top_Xtr: {entry['Top_Xtr']}, Bot_Xtr: {entry['Bot_Xtr']}\n")
        else:
            self.data_display.insert(tk.END, "No data available.\n")
        self.data_display.insert(tk.END, "\n")

    def get_input_values(self):
        return {
            'angle_start': self.angle_start_entry.get(),
            'angle_end': self.angle_end_entry.get(),
            'tip_file': self.tip_file_entry.get(),
            'root_file': self.root_file_entry.get()
        }

    def set_input_values(self, values):
        self.angle_start_entry.delete(0, tk.END)
        self.angle_start_entry.insert(0, values.get('angle_start', self.angle_start_entry.get()))

        self.angle_end_entry.delete(0, tk.END)
        self.angle_end_entry.insert(0, values.get('angle_end', self.angle_end_entry.get()))

        self.tip_file_entry.delete(0, tk.END)
        self.tip_file_entry.insert(0, values.get('tip_file', self.tip_file_entry.get()))

        self.root_file_entry.delete(0, tk.END)
        self.root_file_entry.insert(0, values.get('root_file', self.root_file_entry.get()))



