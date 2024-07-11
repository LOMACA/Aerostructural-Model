from scipy.interpolate import interp1d
from scipy.integrate import quad
import numpy as np
import tkinter as tk
from tkinter import ttk
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from tkinter import filedialog
from airfoil_geometry import AirfoilGeometry


class AirfoilDisplay(ttk.Frame):
    def __init__(self, parent, ltt_display, *args, **kwargs):
        super().__init__(parent, *args, **kwargs)

        self.ltt_display = ltt_display

        self.init_ui()

        self.skin_thickness = 0.001
        self.partition1 = 0.30
        self.partition2 = 0.58

        self.root_airfoil = None
        self.tip_airfoil = None

        self.current_airfoil_type = "Root"
        self.create_output_widgets()

    def create_output_widgets(self):
        output_frame = ttk.Frame(self)
        output_frame.grid(row=1, column=0, sticky="ew", padx=5, pady=5)

        output_frame2 = ttk.Frame(self)
        output_frame2.grid(row=2, column=0, sticky="ew", padx=5, pady=5)

        self.circumference_label = ttk.Label(output_frame, text="Circumference: N/A")
        self.circumference_label.grid(row=0, column=0, padx=5, pady=2)

        # self.leading_edge_circumference_label = ttk.Label(output_frame, text="Leading Edge Circumference: N/A")
        # self.leading_edge_circumference_label.grid(row=0, column=1, padx=5, pady=2)

        self.leading_edge_area_label = ttk.Label(output_frame, text="Leading Edge Area: N/A")
        self.leading_edge_area_label.grid(row=0, column=1, padx=5, pady=2)

        self.middle_section_area_label = ttk.Label(output_frame, text="Middle Section Area: N/A")
        self.middle_section_area_label.grid(row=0, column=2, padx=5, pady=2)

        self.trailing_edge_area_label = ttk.Label(output_frame, text="Trailing Edge Area: N/A")
        self.trailing_edge_area_label.grid(row=0, column=3, padx=5, pady=2)

        self.total_area_label = ttk.Label(output_frame, text="Total Area: N/A")
        self.total_area_label.grid(row=0, column=4, padx=5, pady=2)

        self.distance_to_front_spar_label = ttk.Label(output_frame, text="Distance Quarter-Chord to Front Spar: N/A")
        self.distance_to_front_spar_label.grid(row=0, column=5, padx=5, pady=2)

        self.distance_to_rear_spar_label = ttk.Label(output_frame, text="Distance Quarter-Chord to Rear Spar: N/A")
        self.distance_to_rear_spar_label.grid(row=0, column=6, padx=5, pady=2)

        self.spar_heights_label = ttk.Label(output_frame2, text="Spar Heights: N/A")
        self.spar_heights_label.grid(row=0, column=0, padx=5, pady=2)

        self.spar_distances_label = ttk.Label(output_frame2, text="Spar Distances: N/A")
        self.spar_distances_label.grid(row=0, column=1, padx=5, pady=2)

        self.skin_thickness_label = ttk.Label(output_frame2, text="Skin Thickness: N/A")
        self.skin_thickness_label.grid(row=0, column=2, padx=5, pady=2)

        self.spar_thickness_label = ttk.Label(output_frame2, text="Spar Thickness: N/A")
        self.spar_thickness_label.grid(row=0, column=3, padx=5, pady=2)

        self.max_thickness_label = ttk.Label(output_frame2, text="Max Thickness: N/A")
        self.max_thickness_label.grid(row=0, column=4, padx=5)

    def init_ui(self):

        input_frame = ttk.Frame(self)
        input_frame.grid(row=0, column=0, sticky="ew", padx=5, pady=5)

        self.plot_frame = tk.Frame(self)
        self.plot_frame.grid(row=3, column=0, sticky="nsew", padx=5, pady=5)
        self.plot_frame.grid_columnconfigure(0, weight=1)

        self.fig, self.ax = plt.subplots()
        self.canvas = FigureCanvasTkAgg(self.fig, self.plot_frame)
        self.canvas.get_tk_widget().grid(row=4, column=0, sticky="nsew")

        ttk.Label(input_frame, text="Skin Thickness (m):").grid(row=0, column=0, padx=5)
        self.skin_thickness_entry = ttk.Entry(input_frame, width=5)
        self.skin_thickness_entry.insert(0, "0.001")  # Default to 1 mm
        self.skin_thickness_entry.grid(row=0, column=1, padx=5)

        ttk.Label(input_frame, text="Partition 1 (%):").grid(row=0, column=2, padx=5)
        self.partition1_entry = ttk.Entry(input_frame, width=5)
        self.partition1_entry.insert(0, "30")  # Default 25%
        self.partition1_entry.grid(row=0, column=3, padx=5)

        ttk.Label(input_frame, text="Partition 2 (%):").grid(row=0, column=4, padx=5)
        self.partition2_entry = ttk.Entry(input_frame, width=5)
        self.partition2_entry.insert(0, "58")  # Default 75%
        self.partition2_entry.grid(row=0, column=5, padx=5)

        ttk.Label(input_frame, text="Spar Thickness (m):").grid(row=0, column=6, padx=5)
        self.spar_thickness_entry = ttk.Entry(input_frame, width=5)
        self.spar_thickness_entry.insert(0, "0.007")  # Default 7mm
        self.spar_thickness_entry.grid(row=0, column=7, padx=5)

        self.load_tip_button = ttk.Button(input_frame, text="Load Tip Airfoil", command=self.load_tip_airfoil)
        self.load_tip_button.grid(row=0, column=9, padx=5)

        self.load_root_button = ttk.Button(input_frame, text="Load Root Airfoil", command=self.load_root_airfoil)
        self.load_root_button.grid(row=0, column=10, padx=5)

        self.switch_button = ttk.Button(input_frame, text="Switch Airfoil", command=self.switch_airfoil)
        self.switch_button.grid(row=0, column=13, padx=5)

        self.apst_button = ttk.Button(input_frame, text="Apply Skin Thickness", command=self.apply_skin_thickness)
        self.apst_button.grid(row=0, column=14, padx=5)

        ttk.Button(input_frame, text="Reset Airfoil", command=self.reset_airfoil).grid(row=0, column=18, padx=5)

        self.apply_partitions_button = ttk.Button(input_frame, text="Apply Partitions", command=self.apply_partitions)
        self.apply_partitions_button.grid(row=0, column=16, padx=5)

    def load_root_airfoil(self):
        # Load and rescale the root airfoil
        file_path = filedialog.askopenfilename(title="Select Root Airfoil File", filetypes=[("DAT files", "*.dat")])
        if not file_path:
            return

        self.root_airfoil = AirfoilGeometry(file_path)
        self.root_airfoil.scale_airfoil(self.ltt_display.c_root)  # Scale to root chord length

        if self.current_airfoil_type == "Root":
            self.plot_airfoil("Root Airfoil", self.root_airfoil)

    def load_tip_airfoil(self):

        file_path = filedialog.askopenfilename(title="Select Tip Airfoil File", filetypes=[("DAT files", "*.dat")])
        if not file_path:
            return

        self.tip_airfoil = AirfoilGeometry(file_path)
        self.tip_airfoil.scale_airfoil(self.ltt_display.c_tip)

        if self.current_airfoil_type == "Tip":
            self.plot_airfoil("Tip Airfoil", self.tip_airfoil)

    def apply_skin_thickness(self):

        try:
            self.skin_thickness = float(self.skin_thickness_entry.get())
        except ValueError:
            self.skin_thickness = 0.001  # Default to 1mm if invalid input

        if self.root_airfoil:
            self.root_airfoil.skin_thickness = self.skin_thickness
            self.root_airfoil.apply_skin_thickness()

        if self.tip_airfoil:
            self.tip_airfoil.skin_thickness = self.skin_thickness
            self.tip_airfoil.apply_skin_thickness()  # Apply skin thickness

        if self.current_airfoil_type == "Root" and self.root_airfoil:
            self.plot_airfoil("Root Airfoil", self.root_airfoil)
        elif self.current_airfoil_type == "Tip" and self.tip_airfoil:
            self.plot_airfoil("Tip Airfoil", self.tip_airfoil)

        self.update_coordinates()
        self.update_output_labels()

    def reset_airfoil(self):
        if self.root_airfoil:
            self.root_airfoil.reset_airfoil()
        if self.tip_airfoil:
            self.tip_airfoil.reset_airfoil()

        if self.current_airfoil_type == "Root" and self.root_airfoil:
            self.plot_airfoil("Root Airfoil", self.root_airfoil)
        elif self.current_airfoil_type == "Tip" and self.tip_airfoil:
            self.plot_airfoil("Tip Airfoil", self.tip_airfoil)
        self.update_output_labels()

    def adjust_spars_to_original_coordinates(self, airfoil):

        spar_thickness = float(self.spar_thickness_entry.get())

        for i in range(len(airfoil.coordinates)):
            x, y = airfoil.original_coordinates[i]
            airfoil.coordinates[i] = (x, y + (1 if y >= 0 else -1) * spar_thickness / 2)
        self.update_output_labels()

    def apply_partitions(self):
        try:
            self.partition1 = float(self.partition1_entry.get()) / 100
            self.partition2 = float(self.partition2_entry.get()) / 100
        except ValueError:
            self.partition1 = 0.25
            self.partition2 = 0.6


        if self.root_airfoil:
            self.root_airfoil.update_partitions(self.partition1, self.partition2)

        if self.tip_airfoil:
            self.tip_airfoil.update_partitions(self.partition1, self.partition2)


        self.plot_airfoil("Root Airfoil with Spars", self.root_airfoil)
        self.plot_airfoil("Tip Airfoil with Spars", self.tip_airfoil)

        self.update_coordinates()
        if self.current_airfoil_type == "Root":
            self.update_output_labels()
        else:
            self.update_output_labels()

    def update_coordinates(self):
        if self.root_airfoil:
            self.root_airfoil_parameters = self.get_airfoil_parameters("Root")
        if self.tip_airfoil:
            self.tip_airfoil_parameters = self.get_airfoil_parameters("Tip")

    def calculate_circumference(self, coordinates):

        circumference = 0
        for i in range(1, len(coordinates)):
            x1, y1 = coordinates[i - 1]
            x2, y2 = coordinates[i]
            circumference += np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
        circumference += np.sqrt((coordinates[0][0] - coordinates[-1][0]) ** 2 + (coordinates[0][1] - coordinates[-1][1]) ** 2)
        return circumference

    def calculate_spar_heights(self, airfoil):

        partition1_x = airfoil.chord_length * self.partition1
        partition2_x = airfoil.chord_length * self.partition2

        partition1_y_upper = self.find_y_at_x(partition1_x, airfoil.coordinates)
        partition1_y_lower = self.find_y_at_x(partition1_x, airfoil.coordinates[::-1])
        partition2_y_upper = self.find_y_at_x(partition2_x, airfoil.coordinates)
        partition2_y_lower = self.find_y_at_x(partition2_x, airfoil.coordinates[::-1])

        spar1_height = (partition1_y_upper - partition1_y_lower)
        spar2_height = (partition2_y_upper - partition2_y_lower)

        return {"spar1": spar1_height, "spar2": spar2_height}

    def calculate_spar_distances(self, airfoil):

        partition1_x = airfoil.chord_length * self.partition1
        partition2_x = airfoil.chord_length * self.partition2

        distance_to_first_spar = partition1_x
        distance_between_spars = partition2_x - partition1_x
        distance_from_second_spar_to_trailing_edge = airfoil.chord_length - partition2_x

        return {
            "distance_to_first_spar": distance_to_first_spar,
            "distance_between_spars": distance_between_spars,
            "distance_from_second_spar_to_trailing_edge": distance_from_second_spar_to_trailing_edge
        }

    def get_skin_thickness(self):

        try:
            skin_thickness = float(self.skin_thickness_entry.get())
            return skin_thickness
        except ValueError:
            return 0.001

    def get_spar_thickness(self):

        try:
            spar_thickness = float(self.spar_thickness_entry.get())
            return spar_thickness
        except ValueError:
            return 0.01

    def calculate_distances_to_spars(self, airfoil):
        quarter_chord_x = 0.25 * airfoil.chord_length

        front_spar_x = airfoil.chord_length * self.partition1
        rear_spar_x = airfoil.chord_length * self.partition2

        distance_to_front_spar = front_spar_x - quarter_chord_x
        distance_to_rear_spar = rear_spar_x - quarter_chord_x

        return {
            'distance_to_front_spar': distance_to_front_spar,
            'distance_to_rear_spar': distance_to_rear_spar
        }

    def calculate_leading_edge_circumference(self, coordinates, partition1_x):
        # Find the indices of the leading edge and the first spar position
        leading_edge_index = np.argmin([coord[0] for coord in coordinates])
        spar_index = next(i for i, coord in enumerate(coordinates) if coord[0] >= partition1_x)

        if spar_index > leading_edge_index:
            upper_surface = coordinates[leading_edge_index:spar_index + 1]
            lower_surface = coordinates[spar_index:] + coordinates[:leading_edge_index + 1][::-1]
        else:
            upper_surface = coordinates[:spar_index + 1]
            lower_surface = coordinates[leading_edge_index:][::-1]

        leading_edge_circumference_lower = 0
        leading_edge_circumference_upper = 0

        for i in range(1, len(lower_surface)):
            x1, y1 = lower_surface[i - 1]
            x2, y2 = lower_surface[i]
            leading_edge_circumference_lower += np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)

        for i in range(1, len(upper_surface)):
            x1, y1 = upper_surface[i - 1]
            x2, y2 = upper_surface[i]
            leading_edge_circumference_upper += np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)

        return leading_edge_circumference_lower + leading_edge_circumference_upper

    def calculate_middle_section_circumference(self, airfoil):

        partition1_x = airfoil.chord_length * self.partition1
        partition2_x = airfoil.chord_length * self.partition2

        spar1_index = np.argmax([coord[0] >= partition1_x for coord in airfoil.coordinates])
        spar2_index = np.argmax([coord[0] >= partition2_x for coord in airfoil.coordinates])

        middle_section = airfoil.coordinates[spar1_index:spar2_index + 1]

        middle_section_circumference = 0
        for i in range(1, len(middle_section)):
            x1, y1 = middle_section[i - 1]
            x2, y2 = middle_section[i]
            middle_section_circumference += np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)

        return middle_section_circumference

    def calculate_trailing_edge_circumference(self, airfoil):

        partition2_x = airfoil.chord_length * self.partition2

        spar2_index = np.argmax([coord[0] >= partition2_x for coord in airfoil.coordinates])

        trailing_edge_section = airfoil.coordinates[spar2_index:]

        trailing_edge_circumference = 0
        for i in range(1, len(trailing_edge_section)):
            x1, y1 = trailing_edge_section[i - 1]
            x2, y2 = trailing_edge_section[i]
            trailing_edge_circumference += np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)

        return trailing_edge_circumference

    def interpolate_surfaces(self, coordinates):

        upper_surface = [coord for coord in coordinates if coord[1] >= 0]
        lower_surface = [coord for coord in coordinates if coord[1] < 0]

        upper_surface = np.unique(upper_surface, axis=0)
        lower_surface = np.unique(lower_surface, axis=0)

        upper_surface = upper_surface[np.argsort(upper_surface[:, 0])]
        lower_surface = lower_surface[np.argsort(lower_surface[:, 0])]

        upper_x, upper_y = zip(*upper_surface)
        lower_x, lower_y = zip(*lower_surface)

        upper_interp = interp1d(upper_x, upper_y, kind='quadratic', fill_value="extrapolate")
        lower_interp = interp1d(lower_x, lower_y, kind='quadratic', fill_value="extrapolate")

        return upper_interp, lower_interp

    def calculate_total_area(self, coordinates):
        upper_interp, lower_interp = self.interpolate_surfaces(coordinates)

        def integrand(x):
            return upper_interp(x) - lower_interp(x)

        area, _ = quad(integrand, 0, 1)
        return abs(area)

    def calculate_leading_edge_area(self, total_area):
        return total_area * self.partition1

    def calculate_middle_section_area(self, total_area):
        return total_area * (self.partition2 - self.partition1)

    def calculate_trailing_edge_area(self, total_area):
        return total_area * (1 - self.partition2)

    def plot_airfoil(self, title, airfoil):

        self.ax.clear()

        x_orig = [point[0] for point in airfoil.original_coordinates]
        y_orig = [point[1] for point in airfoil.original_coordinates]
        self.ax.plot(x_orig, y_orig, 'b-', label='Original Airfoil')

        x = [point[0] for point in airfoil.coordinates]
        y = [point[1] for point in airfoil.coordinates]
        self.ax.plot(x, y, 'k-', label='Airfoil with Skin Thickness')

        partition1_x = airfoil.chord_length * self.partition1
        partition2_x = airfoil.chord_length * self.partition2

        partition1_y_upper = self.find_y_at_x(partition1_x, airfoil.coordinates)
        partition1_y_lower = self.find_y_at_x(partition1_x, airfoil.coordinates[::-1])
        partition2_y_upper = self.find_y_at_x(partition2_x, airfoil.coordinates)
        partition2_y_lower = self.find_y_at_x(partition2_x, airfoil.coordinates[::-1])


        partition1_y_upper_scaled = self.find_y_at_x(partition1_x, airfoil.coordinates)
        partition1_y_lower_scaled = self.find_y_at_x(partition1_x, airfoil.coordinates[::-1])
        partition2_y_upper_scaled = self.find_y_at_x(partition2_x, airfoil.coordinates)
        partition2_y_lower_scaled = self.find_y_at_x(partition2_x, airfoil.coordinates[::-1])

        spar_thickness = float(self.spar_thickness_entry.get())
        spar_thickness_visual = spar_thickness * 25

        self.ax.plot([partition1_x, partition1_x], [partition1_y_lower, partition1_y_upper_scaled], 'k-', linewidth=spar_thickness_visual,
                     label=f'Partition 1 at {self.partition1 * 100:.1f}%')
        self.ax.plot([partition2_x, partition2_x], [partition2_y_lower, partition2_y_upper_scaled], 'k-', linewidth=spar_thickness_visual,
                     label=f'Partition 2 at {self.partition2 * 100:.1f}%')

        self.max_thickness_label.configure(
            text=f"Maximum Thickness: {airfoil.max_thickness:.2f} at {airfoil.max_thickness_location * 100:.1f}%"
        )

        self.ax.set_title(title)
        self.ax.set_xlabel("Chord")
        self.ax.set_ylabel("Thickness")
        self.ax.axis("equal")
        self.ax.legend()

        self.canvas.draw()

    def find_y_at_x(self, x_pos, coordinates):
        closest_index = np.argmin([abs(point[0] - x_pos) for point in coordinates])
        return coordinates[closest_index][1]

    def switch_airfoil(self):

        if self.current_airfoil_type == "Root":
            if self.tip_airfoil:
                self.plot_airfoil("Tip Airfoil", self.tip_airfoil)
                self.current_airfoil_type = "Tip"
        else:
            if self.root_airfoil:  #
                self.plot_airfoil("Root Airfoil", self.root_airfoil)
                self.current_airfoil_type = "Root"

        self.canvas.draw()
        self.update_output_labels()

    def get_airfoil_parameters(self, airfoil_type):
        if airfoil_type == "Root" and self.root_airfoil:
            total_area = self.calculate_total_area(self.root_airfoil.coordinates)
            spar_distances = self.calculate_distances_to_spars(self.root_airfoil)
            return {
                "original": self.root_airfoil.original_coordinates,
                "thickened": self.root_airfoil.coordinates,
                "circumference": self.calculate_circumference(self.root_airfoil.original_coordinates),
                "spar_heights": self.calculate_spar_heights(self.root_airfoil),
                "spar_distances": self.calculate_spar_distances(self.root_airfoil),
                "skin_thickness": self.get_skin_thickness(),
                "spar_thickness": self.get_spar_thickness(),
                "leading_edge_area": self.calculate_leading_edge_area(total_area),
                "middle_section_area": self.calculate_middle_section_area(total_area),
                "trailing_edge_area": self.calculate_trailing_edge_area(total_area),
                "total_area": total_area,
                "distance_to_front_spar": spar_distances['distance_to_front_spar'],
                "distance_to_rear_spar": spar_distances['distance_to_rear_spar']
            }

        elif airfoil_type == "Tip" and self.tip_airfoil:
            total_area = self.calculate_total_area(self.tip_airfoil.coordinates)
            spar_distances = self.calculate_distances_to_spars(self.tip_airfoil)
            return {
                "original": self.tip_airfoil.original_coordinates,
                "thickened": self.tip_airfoil.coordinates,
                "circumference": self.calculate_circumference(self.tip_airfoil.coordinates),
                "spar_heights": self.calculate_spar_heights(self.tip_airfoil),
                "spar_distances": self.calculate_spar_distances(self.tip_airfoil),
                "skin_thickness": self.get_skin_thickness(),
                "spar_thickness": self.get_spar_thickness(),
                "leading_edge_area": self.calculate_leading_edge_area(total_area),
                "middle_section_area": self.calculate_middle_section_area(total_area),
                "trailing_edge_area": self.calculate_trailing_edge_area(total_area),
                "total_area": total_area,
                "distance_to_front_spar": spar_distances['distance_to_front_spar'],
                "distance_to_rear_spar": spar_distances['distance_to_rear_spar']
            }

    def update_output_labels(self):
        if self.current_airfoil_type == "Root" and self.root_airfoil:
            params = self.get_airfoil_parameters("Root")
        elif self.current_airfoil_type == "Tip" and self.tip_airfoil:
            params = self.get_airfoil_parameters("Tip")
        else:
            return

        self.circumference_label.config(text=f"Circumference: {np.round(params['circumference'], 2)} m")
        self.leading_edge_area_label.config(text=f"Leading Edge Area: {np.round(params['leading_edge_area'], 4)} m^2")
        self.middle_section_area_label.config(text=f"Middle Section Area: {np.round(params['middle_section_area'], 4)} m^2")
        self.trailing_edge_area_label.config(text=f"Trailing Edge Area: {np.round(params['trailing_edge_area'], 4)} m^2")
        self.spar_heights_label.config(text=f"Spar Heights: {self.format_spar_heights(params['spar_heights'])}")
        self.spar_distances_label.config(text=f"Spar Distances: {self.format_spar_distances(params['spar_distances'])}")
        self.skin_thickness_label.config(text=f"Skin Thickness: {np.round(params['skin_thickness'], 4)} m")
        self.spar_thickness_label.config(text=f"Spar Thickness: {np.round(params['spar_thickness'], 4)} m")
        self.total_area_label.config(text=f"Total Area: {np.round(params['total_area'], 2)} m^2")
        self.distance_to_front_spar_label.config(text=f"Distance Quarter-Chord to Front Spar: {np.round(params['distance_to_front_spar'], 4)} m")
        self.distance_to_rear_spar_label.config(text=f"Distance Quarter-Chord to Rear Spar: {np.round(params['distance_to_rear_spar'], 4)} m")

    def format_spar_heights(self, spar_heights):
        return f"Spar 1: {np.round(spar_heights['spar1'], 4)} m, Spar 2: {np.round(spar_heights['spar2'], 4)} m"

    def format_spar_distances(self, spar_distances):
        return f"To Spar 1: {np.round(spar_distances['distance_to_first_spar'], 4)} m, Between Spars: {np.round(spar_distances['distance_between_spars'], 4)} m, To Trailing Edge: {np.round(spar_distances['distance_from_second_spar_to_trailing_edge'], 4)} m"

    def get_input_values(self):
        return {
            'skin_thickness': self.skin_thickness_entry.get(),
            'partition1': self.partition1_entry.get(),
            'partition2': self.partition2_entry.get(),
            'spar_thickness': self.spar_thickness_entry.get()
        }

    def set_input_values(self, values):
        self.skin_thickness_entry.delete(0, tk.END)
        self.skin_thickness_entry.insert(0, values.get('skin_thickness', "0.001"))

        self.partition1_entry.delete(0, tk.END)
        self.partition1_entry.insert(0, values.get('partition1', "30"))

        self.partition2_entry.delete(0, tk.END)
        self.partition2_entry.insert(0, values.get('partition2', "58"))

        self.spar_thickness_entry.delete(0, tk.END)
        self.spar_thickness_entry.insert(0, values.get('spar_thickness', "0.01"))
