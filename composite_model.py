import numpy as np
import tkinter as tk
from tkinter import ttk
from scipy.optimize import fsolve

class CompositeModel(ttk.Frame):
    def __init__(self, parent, isa_display, ltt_display, airfoil_geometry_display, wing_geometry,
                 structural_model_torenbeek_display, *args, **kwargs):
        super().__init__(parent, *args, **kwargs)

        self.isa_display = isa_display
        self.ltt_display = ltt_display
        self.airfoil_geometry_display = airfoil_geometry_display
        self.wing_geometry = wing_geometry
        self.structural_model_torenbeek_display = structural_model_torenbeek_display
        self.load_case = tk.StringVar()
        self.inertia_relief = tk.StringVar()
        self.create_input_widgets()
        self.show_load_case_specific_inputs()
        self.show_inertia_relief_inputs()

    def create_input_widgets(self):
        input_frame = ttk.LabelFrame(self, text="Inputs")
        input_frame.grid(row=0, column=0, padx=5, pady=5, sticky="nsew")

        load_case_frame = ttk.LabelFrame(input_frame, text="Select Load Case")
        load_case_frame.grid(row=0, column=0, columnspan=2, padx=5, pady=5, sticky="ew")

        self.load_case.set("Maneuver")
        ttk.Radiobutton(load_case_frame, text="Maneuver", variable=self.load_case, value="Maneuver",
                        command=self.show_load_case_specific_inputs).grid(row=0, column=0, padx=5, pady=5)
        ttk.Radiobutton(load_case_frame, text="Gust", variable=self.load_case, value="Gust",
                        command=self.show_load_case_specific_inputs).grid(row=0, column=1, padx=5, pady=5)

        ttk.Label(input_frame, text="Design Maximum Take-Off Mass (kg):").grid(row=1, column=0, padx=2, pady=2)
        self.MTOM_entry = ttk.Entry(input_frame)
        self.MTOM_entry.grid(row=1, column=1, padx=2, pady=2)

        self.spanwise_position_label = ttk.Label(input_frame, text="Spanwise Position (section 0 to n):")
        self.spanwise_position_label.grid(row=2, column=0, padx=2, pady=2)
        self.spanwise_position_entry = ttk.Entry(input_frame)
        self.spanwise_position_entry.grid(row=2, column=1, padx=2, pady=2)
        self.update_spanwise_position_label()

        self.fiber_volume_fraction_skin_label = ttk.Label(input_frame, text="Fiber Volume Fraction (Skin) (%):")
        self.fiber_volume_fraction_skin_label.grid(row=3, column=0, padx=2, pady=2)
        self.fiber_volume_fraction_skin_entry = ttk.Entry(input_frame)
        self.fiber_volume_fraction_skin_entry.grid(row=3, column=1, padx=2, pady=2)

        self.fiber_volume_fraction_spars_label = ttk.Label(input_frame, text="Fiber Volume Fraction (Spars) (%):")
        self.fiber_volume_fraction_spars_label.grid(row=4, column=0, padx=2, pady=2)
        self.fiber_volume_fraction_spars_entry = ttk.Entry(input_frame)
        self.fiber_volume_fraction_spars_entry.grid(row=4, column=1, padx=2, pady=2)

        ttk.Label(input_frame, text="Section CG Position (% chord):").grid(row=5, column=0, padx=2, pady=2)
        self.cg_position_entry = ttk.Entry(input_frame)


        self.load_case_specific_frame = ttk.Frame(input_frame)
        self.load_case_specific_frame.grid(row=6, column=0, columnspan=2, padx=5, pady=5, sticky="ew")

        inertia_relief_frame = ttk.LabelFrame(input_frame, text="Inertia Relief")
        inertia_relief_frame.grid(row=7, column=0, columnspan=2, padx=5, pady=5, sticky="ew")

        self.inertia_relief.set("Yes")
        ttk.Radiobutton(inertia_relief_frame, text="No", variable=self.inertia_relief, value="No",
                        command=self.show_inertia_relief_inputs).grid(row=0, column=0, padx=5, pady=5)
        ttk.Radiobutton(inertia_relief_frame, text="Yes", variable=self.inertia_relief, value="Yes",
                        command=self.show_inertia_relief_inputs).grid(row=0, column=1, padx=5, pady=5)

        self.inertia_relief_frame = ttk.Frame(input_frame)
        self.inertia_relief_frame.grid(row=8, column=0, columnspan=2, padx=5, pady=5, sticky="ew")

        self.result_frame = ttk.LabelFrame(self, text="Results")
        self.result_frame.grid(row=1, column=0, padx=5, pady=5, sticky="ew")

        self.result_display_text = tk.Text(self.result_frame, wrap=tk.WORD, width=60, height=15)
        self.result_display_text.grid(row=0, column=0, padx=5, pady=5)

        ttk.Button(self, text="Compute", command=self.compute).grid(row=2, column=0, padx=5, pady=5, sticky="ew")

    def update_spanwise_position_label(self):
        n = len(self.wing_geometry.wing_geometry)
        self.spanwise_position_label.config(text=f"Spanwise Position (section 0 to {n-1}):")

    def show_load_case_specific_inputs(self):
        for widget in self.load_case_specific_frame.winfo_children():
            widget.destroy()

        if self.load_case.get() == "Gust":
            ttk.Label(self.load_case_specific_frame, text="Derived Gust Velocity (m/s):").grid(row=0, column=0, padx=5,
                                                                                               pady=5)
            self.gust_velocity_entry = ttk.Entry(self.load_case_specific_frame)
            self.gust_velocity_entry.grid(row=0, column=1, padx=5, pady=5)

            ttk.Label(self.load_case_specific_frame, text="Equivalent Airspeed (m/s):").grid(row=1, column=0, padx=5,
                                                                                             pady=5)
            self.equivalent_airspeed_entry = ttk.Entry(self.load_case_specific_frame)
            self.equivalent_airspeed_entry.grid(row=1, column=1, padx=5, pady=5)

    def show_inertia_relief_inputs(self):
        for widget in self.inertia_relief_frame.winfo_children():
            widget.destroy()

        if self.inertia_relief.get() == "Yes":
            ttk.Label(self.inertia_relief_frame, text="Inertia Relief Weight (kg):").grid(row=0, column=0, padx=5,
                                                                                          pady=5)
            self.inertia_relief_weight_entry = ttk.Entry(self.inertia_relief_frame)
            self.inertia_relief_weight_entry.grid(row=0, column=1, padx=5, pady=5)

    def compute(self):
        self.result_display_text.delete("1.0", tk.END)
        MTOM = float(self.MTOM_entry.get())
        spanwise_position = int(self.spanwise_position_entry.get())
        result_text = ""

        if self.load_case.get() == "Maneuver":
            self.load_factor = self.calculate_maneuver_load_factor(MTOM)
            result_text += f"Load Factor: {self.load_factor:.2f}\n"
            self.total_lift = self.compute_total_lift(spanwise_position, self.load_factor)
            result_text += f"Total lift up to spanwise section {spanwise_position}: {self.total_lift:.2f} N\n"

        else:
            gust_velocity = float(self.gust_velocity_entry.get())
            equivalent_airspeed = float(self.equivalent_airspeed_entry.get())
            self.load_factor = self.calculate_gust_load_factor(MTOM, gust_velocity, equivalent_airspeed)
            result_text += f"Load Factor: {self.load_factor:.2f}\n"
            self.total_lift = self.compute_total_lift(spanwise_position, self.load_factor)
            result_text += f"Total lift up to spanwise section {spanwise_position}: {self.total_lift:.2f} N\n"

        if self.inertia_relief.get() == "Yes":
            self.inertia_relief_weight = float(self.inertia_relief_weight_entry.get())
            result_text += f"Inertia Relief Weight: {self.inertia_relief_weight} kg\n"
        else:
            self.inertia_relief_weight = 0

        airfoil, _, _ = self.wing_geometry.wing_geometry[spanwise_position]

        self.spar1_height, self.spar2_height = self.calculate_spar_heights(airfoil, spanwise_position)
        self.section_length = np.sqrt((self.ltt_display.y[-1] - self.ltt_display.y[spanwise_position]) ** 2 + (
                airfoil[-1][2] - airfoil[0][2]) ** 2)
        self.skin_area = self.calculate_skin_area(airfoil) * self.section_length

        fiber_volume_fraction_skin = float(self.fiber_volume_fraction_skin_entry.get()) / 100
        fiber_volume_fraction_spars = float(self.fiber_volume_fraction_spars_entry.get()) / 100

        self.material_density_skin = self.calculate_composite_density(fiber_volume_fraction_skin)
        self.material_density_spars = self.calculate_composite_density(fiber_volume_fraction_spars)

        self.spar_thickness = float(self.airfoil_geometry_display.spar_thickness_entry.get())

        self.spar1_length = np.sqrt((self.ltt_display.y[-1] - self.ltt_display.y[spanwise_position]) ** 2 + (
                airfoil[-1][2] - airfoil[0][2]) ** 2)
        self.spar2_length = self.spar1_length

        self.spar1_weight = self.spar1_length * self.spar1_height * self.spar_thickness * self.material_density_spars * 9.8
        self.spar2_weight = self.spar2_length * self.spar2_height * self.spar_thickness * self.material_density_spars * 9.8

        self.skin_circumference = self.calculate_circumference(airfoil)
        self.skin_weight = self.skin_circumference * self.spar1_length * float(
            self.airfoil_geometry_display.skin_thickness_entry.get()) * self.material_density_skin * 9.8

        self.wing_weight = self.spar1_weight + self.spar2_weight + self.skin_weight
        self.total_weight = self.spar1_weight + self.spar2_weight + self.skin_weight + self.inertia_relief_weight * 9.8

        self.youngs_modulus_skin = self.calculate_youngs_modulus(fiber_volume_fraction_skin, 139360)
        self.youngs_modulus_spars = self.calculate_youngs_modulus(fiber_volume_fraction_spars, 139360)

        self.shear_modulus_skin = self.calculate_shear_modulus(fiber_volume_fraction_skin)
        self.shear_modulus_spars = self.calculate_shear_modulus(fiber_volume_fraction_spars)

        self.torsion_moment = self.calculate_torsion_moment(spanwise_position)
        self.torsion_flows = self.calculate_torsion(airfoil, spanwise_position, self.torsion_moment)
        self.n_xs1 = float(self.torsion_flows[0])
        self.n_xs2 = float(self.torsion_flows[1])
        self.n_xs3 = float(self.torsion_flows[2])
        self.n_xs_spar_web_1 = self.n_xs1 - self.n_xs2
        self.n_xs_spar_web_2 = self.n_xs2 - self.n_xs3

        self.Q_z = self.calculate_shear_force(self.total_lift, self.total_weight)

        self.bending_moment = self.calculate_bending_moment(spanwise_position, self.total_lift, self.total_weight)

        I_y = self.calculate_moment_of_inertia(self.spar1_height, self.section_length)
        self.stresses = self.calculate_stress(self.torsion_flows)
        self.strains = self.calculate_strain(self.stresses)

        self.skin_thickness = float(self.airfoil_geometry_display.skin_thickness_entry.get())
        self.stresses_mm2 = [stress / 1e6 for stress in self.stresses]
        self.stress_spar_web1 = self.n_xs_spar_web_1 / self.skin_thickness
        self.stress_spar_web2 = self.n_xs_spar_web_2 / self.skin_thickness
        self.stress_spar_web1_mm = self.stress_spar_web1 / 1e6
        self.stress_spar_web2_mm = self.stress_spar_web2 / 1e6

        self.shear_flow = self.calculate_shear_flow(spanwise_position, self.total_lift)

        result_text += f"Spar 1 Height: {self.spar1_height:.2f} m\n"
        result_text += f"Spar 2 Height: {self.spar2_height:.2f} m\n"
        result_text += f"Spar Length: {self.spar1_length:.2f} m\n"
        result_text += f"Skin Area: {self.skin_area:.2f} m^2\n"
        result_text += f"Skin Circumference: {self.skin_circumference:.2f} m\n"
        result_text += f"Material Density (Skin): {self.material_density_skin:.2f} kg/m^3\n"
        result_text += f"Material Density (Spars): {self.material_density_spars:.2f} kg/m^3\n"
        result_text += f"Spar 1 Weight: {self.spar1_weight:.2f} N\n"
        result_text += f"Spar 2 Weight: {self.spar2_weight:.2f} N\n"
        result_text += f"Skin Weight: {self.skin_weight:.2f} N\n"
        result_text += f"Wing Weight: {self.wing_weight:.2f} N\n"
        result_text += f"Total Weight: {self.total_weight:.2f} N\n"
        result_text += f"Young's Modulus (Skin): {self.youngs_modulus_skin:.2f} N/mm^2\n"
        result_text += f"Young's Modulus (Spars): {self.youngs_modulus_spars:.2f} N/mm^2\n"
        result_text += f"Shear Modulus (Skin): {self.shear_modulus_skin:.2f} N/mm^2\n"
        result_text += f"Shear Modulus (Spars): {self.shear_modulus_spars:.2f} N/mm^2\n"
        result_text += f"Torsion Moment: {self.torsion_moment:.2f} N*m\n"
        result_text += f"Shear Flows due to Torsion (Nose, Spar, Rear): {self.torsion_flows} N/m \n"
        result_text += f"Shear Flows due to Torsion (Sparweb 1 and 2): {self.n_xs_spar_web_1, self.n_xs_spar_web_2} N/m\n"
        result_text += f"Bending Moment: {self.bending_moment:.2f} N*m\n"
        result_text += f"Moment of Inertia (I_y): {I_y:.6f} m^4\n"
        result_text += f"Stresses (Nose, Spar, Rear): {self.stresses_mm2} N/mm^2 \n"
        result_text += f"Stresses (Sparweb 1, Sparweb 2): {self.stress_spar_web1_mm, self.stress_spar_web2_mm} N/mm^2\n"
        result_text += f"Shear Force: {self.Q_z:.2f} N\n"
        result_text += f"Shear Flow Spar: {self.shear_flow:.2f} N/m \n"

        self.result_display_text.insert(tk.END, result_text)

    def calculate_youngs_modulus(self, fiber_volume_fraction, base_modulus):
        return base_modulus * fiber_volume_fraction / 0.60

    def calculate_shear_modulus(self, fiber_volume_fraction):
        G_m = 1019 #MPa
        G_f = 16000 #MPa
        return G_m * (1 + 0.4 * fiber_volume_fraction ** 0.5) / (
                    (1 - fiber_volume_fraction) ** 1.45 + (G_m / G_f) * fiber_volume_fraction)

    def calculate_maneuver_load_factor(self, MTOM):
        MTOM = float(self.MTOM_entry.get())
        MTOW = MTOM * 9.8
        MTOW_lb = MTOW / 4.448
        n_man = 2.1 + (24000 / (MTOW_lb + 10000))
        return n_man

    def calculate_gust_load_factor(self, MTOM, gust_velocity, equivalent_airspeed):
        U_de = float(self.gust_velocity_entry.get())
        V = float(self.equivalent_airspeed_entry.get())
        rho0 = 1.225
        rho = self.isa_display.rho
        S = self.ltt_display.S
        g = 9.8
        W = float(self.MTOM_entry.get()) * g
        a = self.ltt_display.aw
        mgc = self.wing_geometry.calculate_mgc()
        mug = (2 * (W / S)) / (rho * mgc * a * g)
        kg = (0.88 * mug) / (5.3 + mug)
        n_gust = 1 + (kg * rho0 * U_de * V * a) / (2 * (W / S))
        return n_gust

    def compute_total_lift(self, spanwise_position, load_factor):
        L_dist = self.ltt_display.L_dist
        total_lift = np.sum(L_dist[spanwise_position:]) * load_factor
        return total_lift

    def calculate_spar_heights(self, airfoil, spanwise_position):
        partition1_x = self.airfoil_geometry_display.partition1 * self.wing_geometry.get_chord_length(spanwise_position)
        partition2_x = self.airfoil_geometry_display.partition2 * self.wing_geometry.get_chord_length(spanwise_position)

        partition1_z_upper = self.find_z_at_x(partition1_x, airfoil, upper=True)
        partition1_z_lower = self.find_z_at_x(partition1_x, airfoil, upper=False)
        partition2_z_upper = self.find_z_at_x(partition2_x, airfoil, upper=True)
        partition2_z_lower = self.find_z_at_x(partition2_x, airfoil, upper=False)

        spar1_height = abs(partition1_z_upper - partition1_z_lower)
        spar2_height = abs(partition2_z_upper - partition2_z_lower)

        return spar1_height, spar2_height

    def find_z_at_x(self, x_pos, coordinates, upper=True):
        if upper:
            filtered_coords = [point for point in coordinates if point[1] > np.mean([p[1] for p in coordinates])]
        else:
            filtered_coords = [point for point in coordinates if point[1] <= np.mean([p[1] for p in coordinates])]

        if not filtered_coords:
            return 0

        closest_index = np.argmin([abs(point[0] - x_pos) for point in filtered_coords])
        return filtered_coords[closest_index][2]

    def calculate_skin_area(self, coordinates):
        upper_surface = [coord for coord in coordinates if coord[1] > np.mean([p[1] for p in coordinates])]
        lower_surface = [coord for coord in coordinates if coord[1] <= np.mean([p[1] for p in coordinates])]

        def calculate_surface_area(surface):
            x = np.array([coord[0] for coord in surface])
            z = np.array([coord[2] for coord in surface])
            area = np.trapz(z, x)
            return area

        upper_area = calculate_surface_area(upper_surface)
        lower_area = calculate_surface_area(lower_surface)

        return abs(upper_area) + abs(lower_area)

    def calculate_circumference(self, coordinates):
        upper_surface = [coord for coord in coordinates if coord[1] > np.mean([p[1] for p in coordinates])]
        lower_surface = [coord for coord in coordinates if coord[1] <= np.mean([p[1] for p in coordinates])]

        def calculate_surface_circumference(surface):
            x = np.array([coord[0] for coord in surface])
            z = np.array([coord[2] for coord in surface])
            circumference = np.sum(np.sqrt(np.diff(x) ** 2 + np.diff(z) ** 2))
            return circumference

        upper_circumference = calculate_surface_circumference(upper_surface)
        lower_circumference = calculate_surface_circumference(lower_surface)

        return abs(upper_circumference) + abs(lower_circumference)

    def calculate_composite_density(self, fiber_volume_fraction):
        rho_f = 1780
        rho_m = 1240
        return rho_f * fiber_volume_fraction + rho_m * (1 - fiber_volume_fraction)

    def calculate_torsion(self, airfoil, spanwise_position, torsion_moment):
        nose_thickness = float(self.airfoil_geometry_display.skin_thickness_entry.get())
        spar_thickness = float(self.airfoil_geometry_display.spar_thickness_entry.get())

        spar1_height, spar2_height = self.calculate_spar_heights(airfoil, spanwise_position)

        spar_length = self.wing_geometry.get_chord_length(spanwise_position) * (self.airfoil_geometry_display.partition2 - self.airfoil_geometry_display.partition1)

        nose_length = self.calculate_circumference(airfoil) * self.airfoil_geometry_display.partition1
        trailing_edge_length = self.calculate_circumference(airfoil) * (1 - self.airfoil_geometry_display.partition2)

        G_skin = self.calculate_shear_modulus(float(self.fiber_volume_fraction_skin_entry.get()) / 100)
        G_spars = self.calculate_shear_modulus(float(self.fiber_volume_fraction_spars_entry.get()) / 100)

        A_m1 = self.calculate_skin_area(airfoil) * self.airfoil_geometry_display.partition1
        A_m2 = self.calculate_skin_area(airfoil) * (
                self.airfoil_geometry_display.partition2 - self.airfoil_geometry_display.partition1)
        A_m3 = self.calculate_skin_area(airfoil) * (1 - self.airfoil_geometry_display.partition2)

        def equations(shear_flows):
            n_xs1, n_xs2, n_xs3 = shear_flows

            eq1 = 1 / (2 * A_m1) * (n_xs1 / nose_thickness * nose_length + (n_xs1 - n_xs2) * spar1_height / spar_thickness) - \
                  1 / (2 * A_m2) * (
                          n_xs2 / nose_thickness * spar_length + (n_xs2 - n_xs1) * spar1_height / spar_thickness + (n_xs2 - n_xs3) * spar2_height / spar_thickness)
            eq2 = 1 / (2 * A_m2 ) * (
                    n_xs2 / nose_thickness * spar_length + (n_xs2 - n_xs1) * spar1_height / spar_thickness + (n_xs2 - n_xs3) * spar2_height / spar_thickness) - \
                  1 / (2 * A_m3) * (n_xs3 / nose_thickness * trailing_edge_length + (n_xs2 - n_xs3) * spar2_height / spar_thickness)
            eq3 = torsion_moment - 2 * A_m1 * n_xs1 - 2 * A_m2 * n_xs2 - 2 * A_m3 * n_xs3

            return [eq1, eq2, eq3]

        initial_guess = np.array([1.0, 1.0, 1.0])
        shear_flows = fsolve(equations, initial_guess)

        return shear_flows

    def calculate_torsion_moment(self, spanwise_position):
        load_factor = self.calculate_maneuver_load_factor(float(self.MTOM_entry.get()))
        total_lift = self.compute_total_lift(spanwise_position, load_factor)
        distance_to_cg = float(self.cg_position_entry.get()) / 100 * self.wing_geometry.get_chord_length(spanwise_position)
        torsion_moment = total_lift * distance_to_cg

        return torsion_moment

    def calculate_spar_width(self, spanwise_position):

        partition1_x = self.airfoil_geometry_display.partition1 * self.wing_geometry.get_chord_length(spanwise_position)
        partition2_x = self.airfoil_geometry_display.partition2 * self.wing_geometry.get_chord_length(spanwise_position)
        spar_width = partition2_x - partition1_x
        return spar_width

    def calculate_moment_of_inertia(self, spar_height, spar_width):

        I_y = (spar_width * spar_height ** 3) / 12
        return I_y

    def calculate_bending_moment(self, spanwise_position, total_lift, total_weight):
        y_position = self.ltt_display.y[spanwise_position]
        M = abs((total_lift - total_weight)) * y_position
        return M

    def calculate_strain(self, stresses):
        youngs_modulus_spars = self.calculate_youngs_modulus(float(self.fiber_volume_fraction_spars_entry.get()) / 100,
                                                             139360)
        strains = [stress / youngs_modulus_spars for stress in stresses]
        return strains

    def calculate_stress(self, torsion_flows):
        skin_thickness = float(self.airfoil_geometry_display.skin_thickness_entry.get())
        stresses = [torsion_flow / skin_thickness for torsion_flow in torsion_flows]
        return stresses

    def calculate_shear_force(self, total_lift, total_weight):
        Q_z = total_lift - total_weight
        return Q_z

    def calculate_shear_flow(self, spanwise_position, Q_z):

        spar1_height, spar2_height = self.calculate_spar_heights(self.wing_geometry.wing_geometry[spanwise_position][0],
                                                                 spanwise_position)

        spar_length = self.wing_geometry.get_chord_length(spanwise_position) * ((
                    self.airfoil_geometry_display.partition2 - self.airfoil_geometry_display.partition1)/100)

        youngs_modulus_spars = self.calculate_youngs_modulus(float(self.fiber_volume_fraction_spars_entry.get()) / 100,
                                                             139360)

        I_y = self.calculate_moment_of_inertia(spar1_height, spar_length)

        z = spar1_height / 2
        S = spar_length * ((spar1_height+spar2_height)/2) * z

        shear_force = -Q_z * (youngs_modulus_spars * S) / (youngs_modulus_spars * I_y)

        return shear_force

    def get_input_values(self):
        values = {
            'load_case': self.load_case.get(),
            'MTOM': self.MTOM_entry.get(),
            'spanwise_position': self.spanwise_position_entry.get(),
            'inertia_relief': self.inertia_relief.get(),
            'cg_position': self.cg_position_entry.get(),
            'fiber_volume_fraction_skin': self.fiber_volume_fraction_skin_entry.get(),
            'fiber_volume_fraction_spars': self.fiber_volume_fraction_spars_entry.get()
        }
        if self.load_case.get() == "Gust":
            values.update({
                'gust_velocity': self.gust_velocity_entry.get(),
                'equivalent_airspeed': self.equivalent_airspeed_entry.get()
            })
        if self.inertia_relief.get() == "Yes":
            values.update({
                'inertia_relief_weight': self.inertia_relief_weight_entry.get()
            })
        return values

    def set_input_values(self, values):
        self.load_case.set(values.get('load_case', 'Maneuver'))
        self.MTOM_entry.delete(0, tk.END)
        self.MTOM_entry.insert(0, values.get('MTOM', ''))
        self.spanwise_position_entry.delete(0, tk.END)
        self.spanwise_position_entry.insert(0, values.get('spanwise_position', ''))
        self.cg_position_entry.delete(0, tk.END)
        self.cg_position_entry.insert(0, values.get('cg_position', ''))

        if 'gust_velocity' in values:
            self.gust_velocity_entry.delete(0, tk.END)
            self.gust_velocity_entry.insert(0, values.get('gust_velocity', ''))

        if 'equivalent_airspeed' in values:
            self.equivalent_airspeed_entry.delete(0, tk.END)
            self.equivalent_airspeed_entry.insert(0, values.get('equivalent_airspeed', ''))

        if 'inertia_relief_weight' in values:
            self.inertia_relief_weight_entry.delete(0, tk.END)
            self.inertia_relief_weight_entry.insert(0, values.get('inertia_relief_weight', ''))

        self.fiber_volume_fraction_skin_entry.delete(0, tk.END)
        self.fiber_volume_fraction_skin_entry.insert(0, values.get('fiber_volume_fraction_skin', ''))

        self.fiber_volume_fraction_spars_entry.delete(0, tk.END)
        self.fiber_volume_fraction_spars_entry.insert(0, values.get('fiber_volume_fraction_spars', ''))

        self.show_load_case_specific_inputs()
        self.show_inertia_relief_inputs()

    def get_results(self):
        results = {
            "Parameter": [
                "Load Factor", "Total Lift", "Inertia Relief Weight", "Spar 1 Height", "Spar 2 Height",
                "Skin Area", "Skin Circumference", "Material Density (Skin)", "Material Density (Spars)",
                "Spar 1 Weight", "Spar 2 Weight", "Skin Weight", "Wing Weight", "Total Weight",
                "Young's Modulus (Skin)", "Young's Modulus (Spars)", "Shear Modulus (Skin)", "Shear Modulus (Spars)",
                "Torsion Moment", "Torsion Flows", "Bending Moment", "Strain", "Stresses (Top, Bottom)",
                "Shear Force", "Shear FLow"
            ],
            "Value": [
                self.load_factor, self.total_lift, self.inertia_relief_weight, self.spar1_height, self.spar2_height,
                self.skin_area, self.skin_circumference, self.material_density_skin,
                self.material_density_spars,
                self.spar1_weight, self.spar2_weight, self.skin_weight, self.wing_weight, self.total_weight,
                self.youngs_modulus_skin, self.youngs_modulus_spars, self.shear_modulus_skin, self.shear_modulus_spars,
                self.torsion_moment, self.torsion_flows, self.bending_moment, self.strains, self.stresses,
                self.Q_z, self.shear_flow
            ]
        }
        return results
