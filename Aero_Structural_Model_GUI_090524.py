import tkinter as tk
from tkinter import ttk
from tkinter import filedialog
from isa_display import ISADisplay
from ltt_display import LTTDisplay
from airfoil_aerodynamic_data_display import AirfoilDataDisplay
from structural_model_torenbeek_display import StructuralModelDisplay
from airfoil_geometry_display import AirfoilDisplay
from wing_geometry import WingGeometry
import json
from composite_model import CompositeModel
import pandas as pd

class MainApplication(tk.Tk):
    def __init__(self):
        super().__init__()

        self.title("Aero-structural Model")
        self.geometry("1000x800")

        style = ttk.Style(self)
        style.configure('TNotebook.Tab', font=('Arial', 14))

        self.notebook = ttk.Notebook(self, style='TNotebook')
        self.notebook.pack(fill=tk.BOTH, expand=True)

        home_page = HomePage(self.notebook, self)
        self.notebook.add(home_page, text="Home Page")

        self.isa_display = ISADisplay(self.notebook)
        self.notebook.add(self.isa_display, text="ISA Model")

        self.airfoil_aerodynamic_data_display = AirfoilDataDisplay(self.notebook)
        self.notebook.add(self.airfoil_aerodynamic_data_display, text="Airfoil Aerodynamic Data")

        self.ltt_display = LTTDisplay(self.notebook, self.isa_display, self.airfoil_aerodynamic_data_display)
        self.notebook.add(self.ltt_display, text="Lifting Line Theory")

        self.structural_display = StructuralModelDisplay(self.notebook, self.isa_display, self.ltt_display)
        self.notebook.add(self.structural_display, text="Aluminium Structural Model")

        self.airfoil_geometry_display = AirfoilDisplay(self.notebook, self.ltt_display)
        self.notebook.add(self.airfoil_geometry_display, text="Airfoil Geometry Construction")

        self.wing_geometry = WingGeometry(self.notebook, self.ltt_display, self.airfoil_geometry_display)
        self.notebook.add(self.wing_geometry, text="Wing Geometry Construction")

        self.composite_model = CompositeModel(self.notebook, self.isa_display, self.ltt_display, self.airfoil_geometry_display, self.wing_geometry, self.structural_display)
        self.notebook.add(self.composite_model, text="Composite Structural Model")

    def get_all_input_fields(self):
        input_fields = {
            'ISA': self.isa_display.get_input_values(),
            'AirfoilData': self.airfoil_aerodynamic_data_display.get_input_values(),
            'LTT': self.ltt_display.get_input_values(),
            'Structural': self.structural_display.get_input_values(),
            'AirfoilGeometry': self.airfoil_geometry_display.get_input_values(),
            'WingGeometry': self.wing_geometry.get_input_values(),
            'CompositeModel': self.composite_model.get_input_values()
        }
        return input_fields

    def set_all_input_fields(self, input_fields):
        self.isa_display.set_input_values(input_fields.get('ISA', {}))
        self.airfoil_aerodynamic_data_display.set_input_values(input_fields.get('AirfoilData', {}))
        self.ltt_display.set_input_values(input_fields.get('LTT', {}))
        self.structural_display.set_input_values(input_fields.get('Structural', {}))
        self.airfoil_geometry_display.set_input_values(input_fields.get('AirfoilGeometry', {}))
        self.wing_geometry.set_input_values(input_fields.get('WingGeometry', {}))
        self.composite_model.set_input_values(input_fields.get('CompositeModel', {}))

    def get_all_results(self):
        results = {
            'Structural Aluminium Model': self.structural_display.get_results(),
            'Aerodynamic Model ': self.ltt_display.get_results(),
            'Composite Model': self.composite_model.get_results()
        }
        return results

class HomePage(ttk.Frame):
    def __init__(self, parent, main_app, *args, **kwargs):
        super().__init__(parent, *args, **kwargs)
        self.main_app = main_app

        self.pack(fill=tk.BOTH, expand=True)

        intro_text = (
            "This is an aero-structural model that estimates loads on wings and computes their weight.\n\n"
            "It can be used in the preliminary design phase of aircraft wings.\n\n"
            "The aerodynamic model is based on the description of the Lifting-Line-Theory by Bertin \n\n"
            "and Cummings (2021). The first structural model is based on the description by Torenbeek (2013).\n\n"
            "The second structural model for composite structures is based on the description by Schürmann (2007).\n\n "
            "To begin, you can load an existing study or start a new one.\n\n"
        )
        label = ttk.Label(self, text=intro_text, font=("Arial", 14))
        label.grid(row=0, column=0, padx=20, pady=20, sticky="nsew")

        save_button = ttk.Button(self, text="Save Study", command=self.save_input_fields)
        save_button.grid(row=1, column=0, padx=5, pady=5, sticky="ew")

        load_button = ttk.Button(self, text="Load Study", command=self.load_input_fields)
        load_button.grid(row=2, column=0, padx=5, pady=5, sticky="ew")

        export_button = ttk.Button(self, text="Export Results", command=self.export_results_to_excel)
        export_button.grid(row=3, column=0, padx=5, pady=5, sticky="ew")

    def save_input_fields(self):
        file_path = filedialog.asksaveasfilename(defaultextension=".json", filetypes=[("JSON files", "*.json")])
        if file_path:
            input_fields = self.main_app.get_all_input_fields()
            with open(file_path, 'w') as file:
                json.dump(input_fields, file)

    def load_input_fields(self):
        file_path = filedialog.askopenfilename(defaultextension=".json", filetypes=[("JSON files", "*.json")])
        if file_path:
            try:
                with open(file_path, 'r') as file:
                    input_fields = json.load(file)
                self.main_app.set_all_input_fields(input_fields)
            except FileNotFoundError:
                print("No saved input fields found.")

    def export_results_to_excel(self):
        file_path = filedialog.asksaveasfilename(defaultextension=".xlsx", filetypes=[("Excel files", "*.xlsx")])
        if file_path:
            results = self.main_app.get_all_results()
            with pd.ExcelWriter(file_path) as writer:
                for module_name, result in results.items():
                    df = pd.DataFrame(result)
                    df.to_excel(writer, sheet_name=module_name, index=False)

if __name__ == "__main__":
    app = MainApplication()
    app.mainloop()
