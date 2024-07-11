import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import tkinter as tk
from tkinter import ttk

class WingGeometry(ttk.Frame):
    def __init__(self, parent, ltt_display, airfoil_geometry_display, *args, **kwargs):
        super().__init__(parent, *args, **kwargs)

        self.ltt_display = ltt_display
        self.airfoil_geometry_display = airfoil_geometry_display

        self.create_input_widgets()
        self.wing_geometry = []

    def create_input_widgets(self):
        input_frame = ttk.LabelFrame(self, text="Inputs")
        input_frame.grid(row=0, column=0, padx=5, pady=5, sticky="nsew")

        ttk.Label(input_frame, text="Wing Dihedral Angle (°):").grid(row=0, column=0, padx=2, pady=2)
        self.dihedral_entry = ttk.Entry(input_frame)
        self.dihedral_entry.grid(row=0, column=1, padx=2, pady=2)

        ttk.Label(input_frame, text="Wing Sweep Angles (°) [comma-separated]:").grid(row=1, column=0, padx=2, pady=2)
        self.sweep_entry = ttk.Entry(input_frame)
        self.sweep_entry.grid(row=1, column=1, padx=2, pady=2)

        ttk.Label(input_frame, text="Spanwise Positions for Sweep Angles (m) [comma-separated]:").grid(row=2, column=0, padx=2, pady=2)
        self.span_positions_entry = ttk.Entry(input_frame)
        self.span_positions_entry.grid(row=2, column=1, padx=2, pady=2)

        ttk.Button(input_frame, text="Generate Wing Box Geometry", command=self.generate_wing_geometry).grid(row=3, column=0, columnspan=2, padx=5, pady=5)
        ttk.Button(input_frame, text="Generate Wing Geometry", command=self.create_continuous_skin).grid(row=4, column=0, columnspan=2, padx=5, pady=5)

        self.plot_frame_geometry = ttk.Frame(self)
        self.plot_frame_geometry.grid(row=1, column=0, padx=5, pady=5, sticky="nsew")

        self.plot_frame_skin = ttk.Frame(self)
        self.plot_frame_skin.grid(row=1, column=1, padx=5, pady=5, sticky="nsew")

    def clear_plot_frame(self, frame):
        for widget in frame.winfo_children():
            widget.destroy()

    def parse_sweep_angles(self):
        sweep_angles = list(map(float, self.sweep_entry.get().split(',')))
        span_positions = list(map(float, self.span_positions_entry.get().split(',')))
        return sweep_angles, span_positions

    def generate_wing_geometry(self):
        self.clear_plot_frame(self.plot_frame_geometry)

        wing_geometry = []

        y_positions = self.ltt_display.y

        for y in y_positions:
            airfoil, partition1, partition2 = self.interpolate_airfoil_shape(y)
            if len(airfoil) == 0:
                print(f"Warning: Empty airfoil at y={y}")
            transformed_airfoil = self.apply_transformations(airfoil, y)
            transformed_partition1 = self.apply_transformations([partition1], y)
            transformed_partition2 = self.apply_transformations([partition2], y)
            wing_geometry.append((transformed_airfoil, transformed_partition1, transformed_partition2))

        self.wing_geometry = wing_geometry
        self.visualize_wing_geometry()

    def interpolate_airfoil_shape(self, y):
        root_params = self.airfoil_geometry_display.get_airfoil_parameters("Root")
        tip_params = self.airfoil_geometry_display.get_airfoil_parameters("Tip")

        root_airfoil = root_params["thickened"]
        tip_airfoil = tip_params["thickened"]

        root_partition1 = root_params["spar_distances"]["distance_to_first_spar"]
        root_partition2 = root_params["spar_distances"]["distance_between_spars"] + root_partition1

        tip_partition1 = tip_params["spar_distances"]["distance_to_first_spar"]
        tip_partition2 = tip_params["spar_distances"]["distance_between_spars"] + tip_partition1

        root_chord = self.ltt_display.c_root
        tip_chord = self.ltt_display.c_tip
        span = self.ltt_display.b

        normalized_root_airfoil = [(x / root_chord, y / root_chord) for x, y in root_airfoil]
        normalized_tip_airfoil = [(x / tip_chord, y / tip_chord) for x, y in tip_airfoil]

        chord_length = root_chord - (root_chord - tip_chord) * (y / (span / 2))

        interpolated_airfoil = []
        for root_coord, tip_coord in zip(normalized_root_airfoil, normalized_tip_airfoil):
            x = root_coord[0] + (tip_coord[0] - root_coord[0]) * (y / (span / 2))
            z = root_coord[1] + (tip_coord[1] - root_coord[1]) * (y / (span / 2))
            interpolated_airfoil.append((x * chord_length, z * chord_length))

        partition1 = (root_partition1 + (tip_partition1 - root_partition1) * (y / (span / 2))) * chord_length
        partition2 = (root_partition2 + (tip_partition2 - root_partition2) * (y / (span / 2))) * chord_length

        return interpolated_airfoil, (partition1, 0), (partition2, 0)

    def apply_transformations(self, coordinates, y):
        transformed_coordinates = []

        sweep_angles, span_positions = self.parse_sweep_angles()
        sweep_angle_deg = self.get_sweep_angle(y, sweep_angles, span_positions)

        try:
            dihedral_angle_deg = float(self.dihedral_entry.get())
        except ValueError:
            print("Invalid dihedral angle input. Please enter a valid number.")
            return []

        sweep_angle = np.radians(sweep_angle_deg)
        dihedral_angle = np.radians(dihedral_angle_deg)

        for coord in coordinates:
            x, z = coord
            x_sweep = x + y * np.tan(sweep_angle)
            z_dihedral = z + y * np.tan(dihedral_angle)
            transformed_coordinates.append((x_sweep, y, z_dihedral))

        return transformed_coordinates

    def get_sweep_angle(self, y, sweep_angles, span_positions):
        for i, span_position in enumerate(span_positions):
            if y < span_position:
                return sweep_angles[i]
        return sweep_angles[-1]

    def get_section_area(self, section_index):
        airfoil, _, _ = self.wing_geometry[section_index]
        area = self.calculate_area(airfoil)
        return area

    def calculate_area(self, airfoil):

        if airfoil[0] != airfoil[-1]:
            airfoil.append(airfoil[0])

        x_coords = [coord[0] for coord in airfoil]
        y_coords = [coord[1] for coord in airfoil]

        return 0.5 * np.abs(np.dot(x_coords, np.roll(y_coords, 1)) - np.dot(y_coords, np.roll(x_coords, 1)))

    def get_section_length(self, section_index):
        if section_index < len(self.wing_geometry) - 1:
            _, _, _ = self.wing_geometry[section_index]
            _, _, _ = self.wing_geometry[section_index + 1]
            y1 = self.ltt_display.y[section_index]
            y2 = self.ltt_display.y[section_index + 1]
            return np.abs(y2 - y1)
        else:
            return 0

    def get_chord_length(self, spanwise_position):
        root_chord = self.ltt_display.c_root
        tip_chord = self.ltt_display.c_tip

        chord_length = root_chord + (tip_chord - root_chord) * (spanwise_position / (len(self.ltt_display.y) - 1))
        return chord_length

    def calculate_mgc(self):
        root_chord = self.ltt_display.c_root
        tip_chord = self.ltt_display.c_tip

        mgc = (root_chord + tip_chord) / 2
        return mgc

    def get_spar_height(self, section_index):
        _, partition1, partition2 = self.wing_geometry[section_index]
        return abs(partition1[0][1] - partition2[0][1])

    def get_spar_thickness(self, section_index):
        return float(self.airfoil_geometry_display.spar_thickness_entry.get())

    def visualize_wing_geometry(self):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        skin_thickness = float(self.airfoil_geometry_display.skin_thickness_entry.get())
        spar_thickness_visual = float(self.airfoil_geometry_display.spar_thickness_entry.get()) * 15

        for airfoil, partition1, partition2 in self.wing_geometry:
            airfoil_coords = np.array(airfoil)
            ax.plot(airfoil_coords[:, 0], airfoil_coords[:, 1], airfoil_coords[:, 2])

            if partition1:
                partition1_coords = np.array(partition1)
                print("Partition 1 coordinates:", partition1_coords)
                if partition1_coords.size > 0:
                    spar1_x = [partition1_coords[0, 0], partition1_coords[0, 0]]
                    spar1_y = [partition1_coords[0, 1], partition1_coords[0, 1]]
                    spar1_z = [min(airfoil_coords[:, 2]), max(airfoil_coords[:, 2])]
                    ax.plot(spar1_x, spar1_y, spar1_z, 'k-', linewidth=spar_thickness_visual)

            if partition2:
                partition2_coords = np.array(partition2)
                print("Partition 2 coordinates:", partition2_coords)
                if partition2_coords.size > 0:
                    spar2_x = [partition2_coords[0, 0], partition2_coords[0, 0]]
                    spar2_y = [partition2_coords[0, 1], partition2_coords[0, 1]]
                    spar2_z = [min(airfoil_coords[:, 2]), max(airfoil_coords[:, 2])]
                    ax.plot(spar2_x, spar2_y, spar2_z, 'k-', linewidth=spar_thickness_visual)

        ax.set_xlabel('X axis')
        ax.set_ylabel('Y axis')
        ax.set_zlabel('Z axis')

        z_min, z_max = ax.get_zlim()
        z_range = z_max - z_min
        ax.set_zlim(z_min, z_min + z_range * 2)

        canvas = FigureCanvasTkAgg(fig, master=self.plot_frame_geometry)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    def create_continuous_skin(self):
        self.clear_plot_frame(self.plot_frame_skin)

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        skin_thickness = float(self.airfoil_geometry_display.skin_thickness_entry.get())
        spar_thickness_visual = float(self.airfoil_geometry_display.spar_thickness_entry.get()) * 15

        for section, partition1, partition2 in self.wing_geometry:
            section_coords = np.array(section)
            ax.plot(section_coords[:, 0], section_coords[:, 1], section_coords[:, 2], color='black')
            ax.plot(section_coords[::-1][:, 0], section_coords[::-1][:, 1], section_coords[::-1][:, 2],
                    color='black')

        for i in range(len(self.wing_geometry) - 1):
            sec1, part1_1, part1_2 = self.wing_geometry[i]
            sec2, part2_1, part2_2 = self.wing_geometry[i + 1]

            for j in range(len(sec1) - 1):
                vertices = [
                    sec1[j], sec1[j + 1],
                    sec2[j + 1], sec2[j]
                ]
                poly = Poly3DCollection([vertices], alpha=0.7)
                poly.set_edgecolor('k')
                ax.add_collection3d(poly)

        for airfoil, partition1, partition2 in self.wing_geometry:
            airfoil_coords = np.array(airfoil)

            if partition1:
                partition1_coords = np.array(partition1)
                if partition1_coords.size > 0:
                    spar1_x = [partition1_coords[0, 0], partition1_coords[0, 0]]
                    spar1_y = [partition1_coords[0, 1], partition1_coords[0, 1]]
                    spar1_z = [min(airfoil_coords[:, 2]), max(airfoil_coords[:, 2])]
                    ax.plot(spar1_x, spar1_y, spar1_z, 'k-', linewidth=spar_thickness_visual)

            if partition2:
                partition2_coords = np.array(partition2)
                if partition2_coords.size > 0:
                    spar2_x = [partition2_coords[0, 0], partition2_coords[0, 0]]
                    spar2_y = [partition2_coords[0, 1], partition2_coords[0, 1]]
                    spar2_z = [min(airfoil_coords[:, 2]), max(airfoil_coords[:, 2])]
                    ax.plot(spar2_x, spar2_y, spar2_z, 'k-', linewidth=spar_thickness_visual)


        if self.wing_geometry:
            root_section, _, _ = self.wing_geometry[-1]
            tip_section, _, _ = self.wing_geometry[0]
            root_section_coords = np.array(root_section)
            tip_section_coords = np.array(tip_section)

            for j in range(len(root_section) - 1):
                vertices = [
                    root_section[j], root_section[j + 1],
                    (root_section[j + 1][0], root_section[j + 1][1], 0), (root_section[j][0], root_section[j][1], 0)
                ]
                poly = Poly3DCollection([vertices], alpha=0.7, facecolor='black', edgecolor='k')
                ax.add_collection3d(poly)

            tip_vertices = []
            for j in range(len(tip_section)):
                tip_vertices.append(tip_section[j])

            tip_poly = Poly3DCollection([tip_vertices], alpha=0.7, facecolor='black',
                                        edgecolor='k')
            ax.add_collection3d(tip_poly)

        ax.set_xlabel('X axis')
        ax.set_ylabel('Y axis')
        ax.set_zlabel('Z axis')

        z_min, z_max = ax.get_zlim()
        z_range = z_max - z_min
        ax.set_zlim(z_min, z_min + z_range * 2)

        canvas = FigureCanvasTkAgg(fig, master=self.plot_frame_skin)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    def get_input_values(self):
        return {
            'dihedral_angle': self.dihedral_entry.get(),
            'sweep_angles': self.sweep_entry.get(),
            'span_positions': self.span_positions_entry.get()
        }

    def set_input_values(self, params):
        self.dihedral_entry.delete(0, tk.END)
        self.dihedral_entry.insert(0, params.get('dihedral_angle', ''))

        self.sweep_entry.delete(0, tk.END)
        self.sweep_entry.insert(0, params.get('sweep_angles', ''))

        self.span_positions_entry.delete(0, tk.END)
        self.span_positions_entry.insert(0, params.get('span_positions', ''))

