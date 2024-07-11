import tkinter as tk
from tkinter import ttk
from structural_computation_torenbeek import StructuralComputation

class StructuralModelDisplay(ttk.Frame):
    def __init__(self, parent, isa_display, ltt_display, *args, **kwargs):
        super().__init__(parent, *args, **kwargs)

        self.isa_display = isa_display
        self.ltt_display = ltt_display
        self.create_input_widgets()
        self.create_result_display()

        self.bind_input_clear_events()

    def create_input_widgets(self):

        input_frame = ttk.LabelFrame(self, text="Inputs")
        input_frame.grid(row=0, column=0, padx=5, pady=5, sticky="nsew")

        input_frame2 = ttk.LabelFrame(self, text="Inputs")
        input_frame2.grid(row=0, column=2, padx=5, pady=5, sticky="nsew")

        ttk.Label(input_frame, text="Design Maximum Take-Off Weight (kg):").grid(row=0, column=0, padx=2, pady=2)
        self.MTOM_entry = ttk.Entry(input_frame)
        self.MTOM_entry.grid(row=0, column=1, padx=2, pady=2)

        ttk.Label(input_frame, text="Gross Weight (kg):").grid(row=1, column=0, padx=2, pady=2)
        self.WG_entry = ttk.Entry(input_frame)
        self.WG_entry.grid(row=1, column=1, padx=2, pady=2)

        ttk.Label(input_frame, text="Derived Gust Velocity (m/s):").grid(row=2, column=0, padx=2, pady=2)
        self.U_de_entry = ttk.Entry(input_frame)
        self.U_de_entry.grid(row=2, column=1, padx=2, pady=2)

        ttk.Label(input_frame, text="Mean geometric chord (m):").grid(row=3, column=0, padx=2, pady=2)
        self.MGC_entry = ttk.Entry(input_frame)
        self.MGC_entry.grid(row=3, column=1, padx=2, pady=2)

        ttk.Label(input_frame, text="Smeared Skin Thickness:").grid(row=4, column=0, padx=2, pady=2)
        self.t_sk_entry = ttk.Entry(input_frame)
        self.t_sk_entry.grid(row=4, column=1, padx=2, pady=2)

        ttk.Label(input_frame, text="Section Thickness at Front Spar (m):").grid(row=5, column=0, padx=2, pady=2)
        self.h_fs_entry = ttk.Entry(input_frame)
        self.h_fs_entry.grid(row=5, column=1, padx=2, pady=2)

        ttk.Label(input_frame, text="Section Thickness at Aft Spar (m):").grid(row=6, column=0, padx=2, pady=2)
        self.h_as_entry = ttk.Entry(input_frame)
        self.h_as_entry.grid(row=6, column=1, padx=2, pady=2)

        ttk.Label(input_frame, text="Chord Length of Airfoil (m):").grid(row=7, column=0, padx=2, pady=2)
        self.chord_entry = ttk.Entry(input_frame)
        self.chord_entry.grid(row=7, column=1, padx=2, pady=2)

        ttk.Label(input_frame, text="Maximum Thickness of Airfoil (m):").grid(row=8, column=0, padx=2, pady=2)
        self.t_entry = ttk.Entry(input_frame)
        self.t_entry.grid(row=8, column=1, padx=2, pady=2)

        ttk.Label(input_frame, text="Lateral Coordinate of Wing Group Weight (m):").grid(row=9, column=0, padx=2, pady=2)
        self.y_wg_entry = ttk.Entry(input_frame)
        self.y_wg_entry.grid(row=9, column=1, padx=2, pady=2)

        ttk.Label(input_frame, text="Engine Weight (kg):").grid(row=10, column=0, padx=2, pady=2)
        self.We_entry = ttk.Entry(input_frame)
        self.We_entry.grid(row=10, column=1, padx=2, pady=2)

        ttk.Label(input_frame, text="Lateral Coordinate of Engine Installation (m):").grid(row=11, column=0, padx=2, pady=2)
        self.y_eng_entry = ttk.Entry(input_frame)
        self.y_eng_entry.grid(row=11, column=1, padx=2, pady=2)

        ttk.Label(input_frame, text="Zero Fuel Weight (kg):").grid(row=12, column=0, padx=2, pady=2)
        self.W_zf_entry = ttk.Entry(input_frame)
        self.W_zf_entry.grid(row=12, column=1, padx=2, pady=2)

        ttk.Label(input_frame, text="Length of Fuel Tank (m):").grid(row=13, column=0, padx=2, pady=2)
        self.b_ft_entry = ttk.Entry(input_frame)
        self.b_ft_entry.grid(row=13, column=1, padx=2, pady=2)

        ttk.Label(input_frame, text="Taper Ratio of Fuel Tank:").grid(row=14, column=0, padx=2, pady=2)
        self.TR_ft_entry = ttk.Entry(input_frame)
        self.TR_ft_entry.grid(row=14, column=1, padx=2, pady=2)

        ttk.Label(input_frame, text="Number of Engines:").grid(row=15, column=0, padx=2, pady=2)
        self.N_eng_entry = ttk.Entry(input_frame)
        self.N_eng_entry.grid(row=15, column=1, padx=2, pady=2)

        ttk.Label(input_frame, text="Maximum Landing Weight (kg):").grid(row=16, column=0, padx=2, pady=2)
        self.W_ML_entry = ttk.Entry(input_frame)
        self.W_ML_entry.grid(row=16, column=1, padx=2, pady=2)

        ttk.Label(input_frame, text="Landing Load Factor:").grid(row=17, column=0, padx=2, pady=2)
        self.n_land_entry = ttk.Entry(input_frame)
        self.n_land_entry.grid(row=17, column=1, padx=2, pady=2)

        ttk.Label(input_frame, text="Landing gear fixed to fuselage (y/n)?").grid(row=18, column=0, padx=2, pady=2)
        self.options = ["Yes", "No"]
        self.choice = tk.StringVar()
        self.dropdown = ttk.Combobox(input_frame, textvariable=self.choice, values=self.options, state="readonly")
        self.dropdown.grid(row=18, column=1, padx=2, pady=2)

        ttk.Label(input_frame, text="Mach Number at Design Dive Speed:").grid(row=19, column=0, padx=2, pady=2)
        self.M_D_entry = ttk.Entry(input_frame)
        self.M_D_entry.grid(row=19, column=1, padx=2, pady=2)

        ttk.Label(input_frame, text="Mid-Chord Sweep Angle (°):").grid(row=20, column=0, padx=2, pady=2)
        self.lambda05_entry = ttk.Entry(input_frame)
        self.lambda05_entry.grid(row=20, column=1, padx=2, pady=2)

        ttk.Label(input_frame2, text="Weight of Propulsion System (kg):").grid(row=0, column=0, padx=2, pady=2)
        self.W_pp_entry = ttk.Entry(input_frame2)
        self.W_pp_entry.grid(row=0, column=1, padx=2, pady=2)

        ttk.Label(input_frame2, text="Leading Edge Sweep Angle (°):").grid(row=1, column=0, padx=2, pady=2)
        self.lambda0_entry = ttk.Entry(input_frame2)
        self.lambda0_entry.grid(row=1, column=1, padx=2, pady=2)

        ttk.Label(input_frame2, text="Leading Edge Slats (y/n)?").grid(row=2, column=0, padx=2, pady=2)
        self.options_slats = ["Yes", "No"]
        self.choice_slats = tk.StringVar()
        self.dropdown_slats = ttk.Combobox(input_frame2, textvariable=self.choice_slats, values=self.options_slats, state="readonly")
        self.dropdown_slats.grid(row=2, column=1, padx=2, pady=2)

        ttk.Label(input_frame2, text="Area of fixed leading edge (m^2):").grid(row=3, column=0, padx=2, pady=2)
        self.S_fle_entry = ttk.Entry(input_frame2)
        self.S_fle_entry.grid(row=3, column=1, padx=2, pady=2)

        ttk.Label(input_frame2, text="Area of Slats (m^2):").grid(row=4, column=0, padx=2, pady=2)
        self.S_slat_entry = ttk.Entry(input_frame2)
        self.S_slat_entry.grid(row=4, column=1, padx=2, pady=2)

        ttk.Label(input_frame2, text="Area of Fixed Trailing Edge (m^2):").grid(row=5, column=0, padx=2, pady=2)
        self.S_fte_entry = ttk.Entry(input_frame2)
        self.S_fte_entry.grid(row=5, column=1, padx=2, pady=2)

        ttk.Label(input_frame2, text="Choose Flap Type:").grid(row=6, column=0, padx=2, pady=2)
        self.options_flaps = ["Single or Double Slotted Fowler Flaps", "Triple Slotted Fowler Flaps"]
        self.choice_flaps = tk.StringVar()
        self.dropdown_flaps = ttk.Combobox(input_frame2, textvariable=self.choice_flaps, values=self.options_flaps, state="readonly")
        self.dropdown_flaps.grid(row=6, column=1, padx=2, pady=2)

        ttk.Label(input_frame2, text="Trailing Edge Flap Area (m^2):").grid(row=7, column=0, padx=2, pady=2)
        self.S_tef_entry = ttk.Entry(input_frame2)
        self.S_tef_entry.grid(row=7, column=1, padx=2, pady=2)

        ttk.Label(input_frame2, text="Choose Aileron Type:").grid(row=8, column=0, padx=2, pady=2)
        self.options_ail = ["Unbalanced Ailerons", "Balanced Ailerons", "Mass-Balanced Ailerons"]
        self.choice_ail = tk.StringVar()
        self.dropdown_ail = ttk.Combobox(input_frame2, textvariable=self.choice_ail, values=self.options_ail, state="readonly")
        self.dropdown_ail.grid(row=8, column=1, padx=2, pady=2)

        ttk.Label(input_frame2, text="Aileron Area (m^2):").grid(row=9, column=0, padx=2, pady=2)
        self.S_ail_entry = ttk.Entry(input_frame2)
        self.S_ail_entry.grid(row=9, column=1, padx=2, pady=2)

        ttk.Label(input_frame2, text="Spoiler Area (m^2):").grid(row=10, column=0, padx=2, pady=2)
        self.S_sp_entry = ttk.Entry(input_frame2)
        self.S_sp_entry.grid(row=10, column=1, padx=2, pady=2)

        ttk.Label(input_frame2, text="Winglet Length (m):").grid(row=11, column=0, padx=2, pady=2)
        self.l_tip_entry = ttk.Entry(input_frame2)
        self.l_tip_entry.grid(row=11, column=1, padx=2, pady=2)

        ttk.Label(input_frame2, text="Winglet Area (m^2):").grid(row=12, column=0, padx=2, pady=2)
        self.S_tip_entry = ttk.Entry(input_frame2)
        self.S_tip_entry.grid(row=12, column=1, padx=2, pady=2)

        ttk.Label(input_frame2, text="Total Number of General Ribs:").grid(row=13, column=0, padx=2, pady=2)
        self.N_rib_entry = ttk.Entry(input_frame2)
        self.N_rib_entry.grid(row=13, column=1, padx=2, pady=2)

        ttk.Label(input_frame2, text="Chord Length of the Wing Box (m):").grid(row=14, column=0, padx=2, pady=2)
        self.c_box_entry = ttk.Entry(input_frame2)
        self.c_box_entry.grid(row=14, column=1, padx=2, pady=2)

        ttk.Label(input_frame2, text="Rib Pitch (m):").grid(row=15, column=0, padx=2, pady=2)
        self.r_entry = ttk.Entry(input_frame2)
        self.r_entry.grid(row=15, column=1, padx=2, pady=2)

        ttk.Label(input_frame2, text="Cant Angle of the Winglet (°):").grid(row=16, column=0, padx=2, pady=2)
        self.gamma_tip_entry = ttk.Entry(input_frame2)
        self.gamma_tip_entry.grid(row=16, column=1, padx=2, pady=2)

        ttk.Label(input_frame2, text="Wingtip Span (m):").grid(row=17, column=0, padx=2, pady=2)
        self.b_tip_entry = ttk.Entry(input_frame2)
        self.b_tip_entry.grid(row=17, column=1, padx=2, pady=2)

        ttk.Label(input_frame2, text="Centre Section Span (m):").grid(row=18, column=0, padx=2, pady=2)
        self.b_cs_entry = ttk.Entry(input_frame2)
        self.b_cs_entry.grid(row=18, column=1, padx=2, pady=2)

        ttk.Label(input_frame2, text="Chord Length at Thickness Kink (m):").grid(row=19, column=0, padx=2, pady=2)
        self.c_tk_entry = ttk.Entry(input_frame2)
        self.c_tk_entry.grid(row=19, column=1, padx=2, pady=2)

        ttk.Label(input_frame2, text="Thickness at Thickness Kink (m):").grid(row=20, column=0, padx=2, pady=2)
        self.t_tk_entry = ttk.Entry(input_frame2)
        self.t_tk_entry.grid(row=20, column=1, padx=2, pady=2)

        ttk.Label(input_frame2, text="Leading Edge made of CFC (y/n)?:").grid(row=21, column=0, padx=2, pady=2)
        self.options_fle = ["Yes", "No"]
        self.choice_fle = tk.StringVar()
        self.dropdown_fle = ttk.Combobox(input_frame2, textvariable=self.choice_fle, values=self.options_fle, state="readonly")
        self.dropdown_fle.grid(row=21, column=1, padx=2, pady=2)

        ttk.Label(input_frame2, text="Trailing Edge Flaps made of CFC (y/n)?:").grid(row=22, column=0, padx=2, pady=2)
        self.options_tef = ["Yes", "No"]
        self.choice_tef = tk.StringVar()
        self.dropdown_tef = ttk.Combobox(input_frame2, textvariable=self.choice_tef, values=self.options_tef, state="readonly")
        self.dropdown_tef.grid(row=22, column=1, padx=2, pady=2)

        ttk.Label(input_frame2, text="Ailerons made of CFC (y/n)?:").grid(row=23, column=0, padx=2, pady=2)
        self.options_ail_cfc = ["Yes", "No"]
        self.choice_ail_cfc = tk.StringVar()
        self.dropdown_ail_cfc = ttk.Combobox(input_frame2, textvariable=self.choice_ail_cfc, values=self.options_ail_cfc, state="readonly")
        self.dropdown_ail_cfc.grid(row=23, column=1, padx=2, pady=2)

        ttk.Label(input_frame2, text="Spoilers made of CFC (y/n)?:").grid(row=24, column=0, padx=2, pady=2)
        self.options_spoil = ["Yes", "No"]
        self.choice_spoil = tk.StringVar()
        self.dropdown_spoil = ttk.Combobox(input_frame2, textvariable=self.choice_spoil, values=self.options_spoil, state="readonly")
        self.dropdown_spoil.grid(row=24, column=1, padx=2, pady=2)

        ttk.Label(input_frame2, text="Wing-Fuselage Fairing made of CFC (y/n)?:").grid(row=25, column=0, padx=2, pady=2)
        self.options_w_fus = ["Yes", "No"]
        self.choice_w_fus = tk.StringVar()
        self.dropdown_w_fus = ttk.Combobox(input_frame2, textvariable=self.choice_w_fus, values=self.options_w_fus, state="readonly")
        self.dropdown_w_fus.grid(row=25, column=1, padx=2, pady=2)

        ttk.Button(input_frame2, text="Compute", command=self.compute_weight).grid(row=26, column=0, columnspan=2, padx=5, pady=5)

    def create_result_display(self):
        result_frame = ttk.LabelFrame(self, text="Results")
        result_frame.grid(row=0, column=3, padx=5, pady=0, sticky="nsew")

        self.result_display_text = tk.Text(result_frame, wrap=tk.WORD, width=60, height=40)
        self.result_display_text.grid(row=1, column=0, padx=1, pady=1)

    def bind_input_clear_events(self):
        input_widgets = [
            self.MTOM_entry,
            self.WG_entry,
            self.U_de_entry,
            self.MGC_entry,
            self.t_sk_entry,
            self.h_fs_entry,
            self.h_as_entry,
            self.chord_entry,
            self.t_entry,
            self.y_wg_entry,
            self.We_entry,
            self.y_eng_entry,
            self.W_zf_entry,
            self.b_ft_entry,
            self.TR_ft_entry,
            self.N_eng_entry,
            self.W_ML_entry,
            self.n_land_entry,
            self.M_D_entry,
            self.lambda05_entry,
            self.lambda0_entry,
            self.S_fle_entry,
            self.S_slat_entry,
            self.S_fte_entry,
            self.S_tef_entry,
            self.S_ail_entry,
            self.S_sp_entry,
            self.l_tip_entry,
            self.S_tip_entry,
            self.N_rib_entry,
            self.c_box_entry,
            self.r_entry,
            self.gamma_tip_entry,
            self.b_tip_entry,
            self.b_cs_entry,
            self.c_tk_entry,
            self.t_tk_entry,
            self.dropdown,
            self.dropdown_slats,
            self.dropdown_flaps,
            self.dropdown_ail,
            self.dropdown_ail_cfc,
            self.dropdown_spoil,
            self.dropdown_w_fus]

        for widget in input_widgets:
            widget.bind("<KeyRelease>", self.clear_results)
            if isinstance(widget, ttk.Combobox):
                widget.bind("<<ComboboxSelected>>", self.clear_results)

    def clear_results(self, event=None):
        self.result_display_text.delete("1.0", tk.END)

    def compute_weight(self):

        MTOM = float(self.MTOM_entry.get())
        W_pp = float(self.W_pp_entry.get())
        MG = float(self.WG_entry.get())
        U_de = float(self.U_de_entry.get())
        t_sk = float(self.t_sk_entry.get())
        h_fs = float(self.h_fs_entry.get())
        h_as = float(self.h_as_entry.get())
        chord = float(self.chord_entry.get())
        t = float(self.t_entry.get())
        y_wg = float(self.y_wg_entry.get())
        We = float(self.We_entry.get())
        y_eng = float(self.y_eng_entry.get())
        W_zf = float(self.W_zf_entry.get())
        b_ft = float(self.b_ft_entry.get())
        TR_ft = float(self.TR_ft_entry.get())
        N_eng = int(self.N_eng_entry.get())
        W_ML = float(self.W_ML_entry.get())
        n_land = float(self.n_land_entry.get())
        M_D = float(self.M_D_entry.get())
        lambda_05 = float(self.lambda05_entry.get())
        lambda_0 = float(self.lambda0_entry.get())
        LE_slat_choice = str(self.dropdown_slats.get())
        S_fle = float(self.S_fle_entry.get())
        S_slat = float(self.S_slat_entry.get())
        S_fte = float(self.S_fte_entry.get())
        TE_flap_choice = self.dropdown_flaps.get()
        S_tef = float(self.S_tef_entry.get())
        ail_choice = self.dropdown_ail.get()
        S_ail = float(self.S_ail_entry.get())
        S_sp = float(self.S_sp_entry.get())
        l_tip = float(self.l_tip_entry.get())
        S_tip = float(self.S_tip_entry.get())
        N_rib = int(self.N_rib_entry.get())
        c_box = float(self.c_box_entry.get())
        r = float(self.r_entry.get())
        gamma_tip = float(self.gamma_tip_entry.get())
        b_tip = float(self.b_tip_entry.get())
        b_cs = float(self.b_cs_entry.get())
        c_tk = float(self.c_tk_entry.get())
        t_tk = float(self.t_tk_entry.get())
        choice = str(self.dropdown.get())
        choice_fle = str(self.dropdown_fle.get())
        choice_tef = str(self.dropdown_tef.get())
        choice_ail_cfc = str(self.dropdown_ail_cfc.get())
        choice_spoil = str(self.dropdown_spoil.get())
        choice_w_fus = str(self.dropdown_w_fus.get())

        AR = self.ltt_display.AR
        TR = self.ltt_display.TR
        L = self.ltt_display.L
        rho = self.isa_display.rho
        S = self.ltt_display.S
        a = self.ltt_display.aw
        b = self.ltt_display.b
        v_inf = self.ltt_display.v_inf
        y = self.ltt_display.y
        c_root = self.ltt_display.c_root

        weight = StructuralComputation(MTOM, MG, rho, S, b, a, U_de, v_inf, y, TR, AR,t_sk, t, h_as, h_fs, y_wg, We, y_eng, W_zf, b_ft,
                TR_ft, N_eng, W_pp, W_ML, n_land, M_D, lambda_05, lambda_0, LE_slat_choice, S_fle, S_slat, S_fte,
                TE_flap_choice, S_tef, ail_choice, S_ail, S_sp, l_tip, S_tip, N_rib, c_box, r, gamma_tip, b_tip, b_cs,
                c_tk, t_tk, choice_fle, choice_tef, choice_ail_cfc, choice_spoil, choice_w_fus, c_root, L, choice, chord)

        results_text = weight.compute_weight()

        self.result_display_text.delete("1.0", tk.END)

        self.result_display_text.insert(tk.END, results_text)

    def get_input_values(self):
        return {
            'MTOM': self.MTOM_entry.get(),
            'WG': self.WG_entry.get(),
            'U_de': self.U_de_entry.get(),
            'MGC': self.MGC_entry.get(),
            't_sk': self.t_sk_entry.get(),
            'h_fs': self.h_fs_entry.get(),
            'h_as': self.h_as_entry.get(),
            'chord': self.chord_entry.get(),
            't': self.t_entry.get(),
            'y_wg': self.y_wg_entry.get(),
            'We': self.We_entry.get(),
            'y_eng': self.y_eng_entry.get(),
            'W_zf': self.W_zf_entry.get(),
            'b_ft': self.b_ft_entry.get(),
            'TR_ft': self.TR_ft_entry.get(),
            'N_eng': self.N_eng_entry.get(),
            'W_ML': self.W_ML_entry.get(),
            'n_land': self.n_land_entry.get(),
            'gear_choice': self.dropdown.get(),
            'M_D': self.M_D_entry.get(),
            'lambda05': self.lambda05_entry.get(),
            'W_pp': self.W_pp_entry.get(),
            'lambda0': self.lambda0_entry.get(),
            'slats_choice': self.dropdown_slats.get(),
            'S_fle': self.S_fle_entry.get(),
            'S_slat': self.S_slat_entry.get(),
            'S_fte': self.S_fte_entry.get(),
            'flaps_choice': self.dropdown_flaps.get(),
            'S_tef': self.S_tef_entry.get(),
            'ail_choice': self.dropdown_ail.get(),
            'S_ail': self.S_ail_entry.get(),
            'S_sp': self.S_sp_entry.get(),
            'l_tip': self.l_tip_entry.get(),
            'S_tip': self.S_tip_entry.get(),
            'N_rib': self.N_rib_entry.get(),
            'c_box': self.c_box_entry.get(),
            'r': self.r_entry.get(),
            'gamma_tip': self.gamma_tip_entry.get(),
            'b_tip': self.b_tip_entry.get(),
            'b_cs': self.b_cs_entry.get(),
            'c_tk': self.c_tk_entry.get(),
            't_tk': self.t_tk_entry.get(),
            'fle_choice': self.dropdown_fle.get(),
            'tef_choice': self.dropdown_tef.get(),
            'ail_cfc_choice': self.dropdown_ail_cfc.get(),
            'spoil_choice': self.dropdown_spoil.get(),
            'w_fus_choice': self.dropdown_w_fus.get(),
        }

    def set_input_values(self, params):
        self.MTOM_entry.delete(0, tk.END)
        self.MTOM_entry.insert(0, params.get('MTOM', ''))

        self.WG_entry.delete(0, tk.END)
        self.WG_entry.insert(0, params.get('WG', ''))

        self.U_de_entry.delete(0, tk.END)
        self.U_de_entry.insert(0, params.get('U_de', ''))

        self.MGC_entry.delete(0, tk.END)
        self.MGC_entry.insert(0, params.get('MGC', ''))

        self.t_sk_entry.delete(0, tk.END)
        self.t_sk_entry.insert(0, params.get('t_sk', ''))

        self.h_fs_entry.delete(0, tk.END)
        self.h_fs_entry.insert(0, params.get('h_fs', ''))

        self.h_as_entry.delete(0, tk.END)
        self.h_as_entry.insert(0, params.get('h_as', ''))

        self.chord_entry.delete(0, tk.END)
        self.chord_entry.insert(0, params.get('chord', ''))

        self.t_entry.delete(0, tk.END)
        self.t_entry.insert(0, params.get('t', ''))

        self.y_wg_entry.delete(0, tk.END)
        self.y_wg_entry.insert(0, params.get('y_wg', ''))

        self.We_entry.delete(0, tk.END)
        self.We_entry.insert(0, params.get('We', ''))

        self.y_eng_entry.delete(0, tk.END)
        self.y_eng_entry.insert(0, params.get('y_eng', ''))

        self.W_zf_entry.delete(0, tk.END)
        self.W_zf_entry.insert(0, params.get('W_zf', ''))

        self.b_ft_entry.delete(0, tk.END)
        self.b_ft_entry.insert(0, params.get('b_ft', ''))

        self.TR_ft_entry.delete(0, tk.END)
        self.TR_ft_entry.insert(0, params.get('TR_ft', ''))

        self.N_eng_entry.delete(0, tk.END)
        self.N_eng_entry.insert(0, params.get('N_eng', ''))

        self.W_ML_entry.delete(0, tk.END)
        self.W_ML_entry.insert(0, params.get('W_ML', ''))

        self.n_land_entry.delete(0, tk.END)
        self.n_land_entry.insert(0, params.get('n_land', ''))

        self.dropdown.set(params.get('gear_choice', ''))

        self.M_D_entry.delete(0, tk.END)
        self.M_D_entry.insert(0, params.get('M_D', ''))

        self.lambda05_entry.delete(0, tk.END)
        self.lambda05_entry.insert(0, params.get('lambda05', ''))

        self.W_pp_entry.delete(0, tk.END)
        self.W_pp_entry.insert(0, params.get('W_pp', ''))

        self.lambda0_entry.delete(0, tk.END)
        self.lambda0_entry.insert(0, params.get('lambda0', ''))

        self.dropdown_slats.set(params.get('slats_choice', ''))

        self.S_fle_entry.delete(0, tk.END)
        self.S_fle_entry.insert(0, params.get('S_fle', ''))

        self.S_slat_entry.delete(0, tk.END)
        self.S_slat_entry.insert(0, params.get('S_slat', ''))

        self.S_fte_entry.delete(0, tk.END)
        self.S_fte_entry.insert(0, params.get('S_fte', ''))

        self.dropdown_flaps.set(params.get('flaps_choice', ''))

        self.S_tef_entry.delete(0, tk.END)
        self.S_tef_entry.insert(0, params.get('S_tef', ''))

        self.dropdown_ail.set(params.get('ail_choice', ''))

        self.S_ail_entry.delete(0, tk.END)
        self.S_ail_entry.insert(0, params.get('S_ail', ''))

        self.S_sp_entry.delete(0, tk.END)
        self.S_sp_entry.insert(0, params.get('S_sp', ''))

        self.l_tip_entry.delete(0, tk.END)
        self.l_tip_entry.insert(0, params.get('l_tip', ''))

        self.S_tip_entry.delete(0, tk.END)
        self.S_tip_entry.insert(0, params.get('S_tip', ''))

        self.N_rib_entry.delete(0, tk.END)
        self.N_rib_entry.insert(0, params.get('N_rib', ''))

        self.c_box_entry.delete(0, tk.END)
        self.c_box_entry.insert(0, params.get('c_box', ''))

        self.r_entry.delete(0, tk.END)
        self.r_entry.insert(0, params.get('r', ''))

        self.gamma_tip_entry.delete(0, tk.END)
        self.gamma_tip_entry.insert(0, params.get('gamma_tip', ''))

        self.b_tip_entry.delete(0, tk.END)
        self.b_tip_entry.insert(0, params.get('b_tip', ''))

        self.b_cs_entry.delete(0, tk.END)
        self.b_cs_entry.insert(0, params.get('b_cs', ''))

        self.c_tk_entry.delete(0, tk.END)
        self.c_tk_entry.insert(0, params.get('c_tk', ''))

        self.t_tk_entry.delete(0, tk.END)
        self.t_tk_entry.insert(0, params.get('t_tk', ''))

        self.dropdown_fle.set(params.get('fle_choice', ''))

        self.dropdown_tef.set(params.get('tef_choice', ''))

        self.dropdown_ail_cfc.set(params.get('ail_cfc_choice', ''))

        self.dropdown_spoil.set(params.get('spoil_choice', ''))

        self.dropdown_w_fus.set(params.get('w_fus_choice', ''))

    def get_results(self):
        results = {
            "Parameter": [
                "Maneuver Load Factor", "Gust Load Factor", "Limit Load Factor", "Ultimate Load Factor",
                "Ideal Weight", "Non-Ideal Weight", "Total Wing Weight"
            ],
            "Value": []
        }

        results_text = self.result_display_text.get("1.0", tk.END).strip().split("\n")
        for line in results_text:
            value = line.split(":")[-1].strip()
            results["Value"].append(value)

        return results
