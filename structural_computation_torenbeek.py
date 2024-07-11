import numpy as np

class StructuralComputation:
    def __init__(self, MTOM, MG, rho, S, b, a, U_de, v_inf, y, TR, AR,t_sk, t, h_as, h_fs, y_wg, We, y_eng, W_zf, b_ft,
                TR_ft, N_eng, W_pp, W_ML, n_land, M_D, lambda_05, lambda_0, LE_slat_choice, S_fle, S_slat, S_fte,
                TE_flap_choice, S_tef, ail_choice, S_ail, S_sp, l_tip, S_tip, N_rib, c_box, r, gamma_tip, b_tip, b_cs,
                c_tk, t_tk, choice_fle, choice_tef, choice_ail_cfc, choice_spoil, choice_w_fus, c_root, L, choice, chord):

        self.MTOM = float(MTOM)
        self.MG = float(MG)
        self.g = 9.81
        self.rho = float(rho)
        self.S = float(S)
        self.b = float(b)
        self.a = float(a)
        self.U_de = float(U_de)
        self.V_EAS = float(v_inf)
        self.y = y
        self.TR = float(TR)
        self.AR = float(AR)
        self.t_sk = float(t_sk)
        self.c0 = float(c_root)
        self.t = float(t)
        self.h_as = float(h_as)
        self.h_fs = float(h_fs)
        self.y_wg = float(y_wg)
        self.We = float(We)
        self.y_eng = float(y_eng)
        self.W_Zf = float(W_zf)
        self.b_ft = float(b_ft)
        self.TR_ft = float(TR_ft)
        self.N_eng = float(N_eng)
        self.W_pp = float(W_pp)
        self.W_ML = float(W_ML)
        self.n_land = float(n_land)
        self.M_D = float(M_D)
        self.lambda_05 = float(lambda_05)
        self.lambda_0 = float(lambda_0)
        self.LE_slat_choice = str(LE_slat_choice)
        self.S_fle = float(S_fle)
        self.S_slat = float(S_slat)
        self.S_fte = float(S_fte)
        self.TE_flap_choice = str(TE_flap_choice)
        self.S_tef = float(S_tef)
        self.ail_choice = str(ail_choice)
        self.S_ail = float(S_ail)
        self.S_sp = float(S_sp)
        self.l_tip = float(l_tip)
        self.S_tip = float(S_tip)
        self.N_rib = float(N_rib)
        self.c_box = float(c_box)
        self.r = float(r)
        self.gamma_tip = float(gamma_tip)
        self.b_tip = float(b_tip)
        self.b_cs = float(b_cs)
        self.c_tk = float(c_tk)
        self.t_tk = float(t_tk)
        self.choice_fle = str(choice_fle)
        self.choice_tef = str(choice_tef)
        self.choice_ail_cfc = str(choice_ail_cfc)
        self.choice_spoil = str(choice_spoil)
        self.choice_w_fus = str(choice_w_fus)
        self.c0 = float(c_root)
        self.L = float(L)
        self.choice = str(choice)
        self.chord = float(chord)

    def compute_weight(self):

        results_text = ""

        rho_g = 2796  # specific weight of aluminium alloy
        I2_t = 0.36  # tension second moment of inertia
        I2_c = 0.45  # compression second moment of inertia
        rho_sl = self.rho
        rho_sl_imp = self.rho * 0.062428
        MTOM_lb = self.MTOM / 0.4536  # MTOW in lb
        WG_lb = self.MG / 0.4536  # gross TO weight in lb
        MTOW = self.MTOM * self.g
        MTOW_lb = self.MTOM * (self.g / 0.3048)
        WG = self.MG * self.g
        WG_lb = self.MG * (self.g / 0.3048)
        b_imp = self.b / 0.3048
        rho_imp = self.rho * 0.062428
        g_imp = self.g * 3.28084
        S_imp = self.S * 10.7639
        U_de_imp = self.U_de * 3.28084
        V_EAS_imp = self.V_EAS * 3.28084
        gamma_tip_rad = self.gamma_tip * (np.pi / 180)
        b_tip_half = self.b_tip / 2

        n_man = 2.1 + (24000 / (
                    MTOW_lb + 10000))  # limit load factor according to CS23 for maneuvers - must not exceed 3.8
        mu_g = (2 * WG_lb * b_imp) / (rho_imp * g_imp * S_imp ** 2) * (self.a) ** -1
        K_g = (0.88 * mu_g) / (5.3 + mu_g)
        n_gust = 1 + K_g * ((0.5 * rho_sl_imp * V_EAS_imp * S_imp) / (WG_lb)) * self.a
        ref = self.a * (0.44 * ((rho_sl_imp * U_de_imp * V_EAS_imp) / (n_man - 1)) - 2.65 * ((rho_imp * g_imp * S_imp) / b_imp))
        wing_loading = MTOW_lb / S_imp

        if n_gust > n_man:
            n_limit = n_gust
        else:
            n_limit = n_man

        n_ult = 1.5 * n_limit

        results_text += f"Maneuver Load Factor: {n_man}\n"
        results_text += f"Gust Load Factor: {n_gust}\n"
        results_text += f"Limit Load Factor: {n_limit}\n"
        results_text += f"Ultimate Load Factor: {n_ult}\n"

        lambda_ea = 0 * (np.pi / 180)  # sweep angle of the elastic axis
        b_st = self.b / np.cos(lambda_ea)
        eta_val = self.y / (self.b * 0.5)
        gamma_tor = (1 - eta_val * (1 - self.TR)) * (2 / (1 + self.TR))
        eta_cp = 1 / (3 * n_ult) * (4 / np.pi + (n_ult - 1) * ((1 + 2 * self.TR) / (1 + self.TR)))

        I1 = (1 - eta_val) ** (3 - 2 * self.TR + self.TR ** 2)

        M_BL = (self.L / 4) * eta_cp * b_st * I1
        M_BL_0 = (self.L / 4) * b_st * eta_cp

        sigma_t = 425 * 10 ** 6
        sigma_c = 420 * 10 ** 6

        R_cant = (self.AR * (1 + self.TR)) / (4 * (self.t / self.c0) * np.cos(lambda_ea))  # cantilever ratio
        sigma_r = 1 / (0.5 * ((1 / sigma_t) + (1.25 / sigma_c)))
        eta_t = np.sqrt((1 + (self.h_fs / self.t) ** 2 + (self.h_as / self.t) ** 2) / 3) - self.t_sk / self.t
        W_BL = I2_t * (rho_g / sigma_r) * n_ult * MTOW * R_cant * (eta_cp / eta_t) * b_st
        W_BL_alt = 0.5 * n_ult * WG * b_st * (eta_cp / eta_t) * (b_st / 2 * self.t) * ((rho_g / sigma_t) * I2_t + (rho_g / sigma_c) * I2_c)

        # sigma_u = 4000  # shear stress in upper cover panels
        # sigma_l = 4000  # shear stress in lower cover panels
        sigma = 0.5 * (sigma_t + sigma_c)
        W_SL = 1.2 * n_ult * WG * eta_cp * b_st * (rho_g / sigma)

        t_r = self.t
        t_t = t_r - 0.05
        t_ref = 1
        k_rib = 0.5 * 10 ** -3
        W_rib = self.N_rib*k_rib * rho_g * self.S * (t_ref + ((t_r + t_t) / 2))
        W_w = W_rib + W_SL + W_BL

        b_half = self.b / 2
        y_cp = eta_cp * b_half

        R_in_wg = W_w / self.MTOM  # relief factor on the root bending moment
        eta_eng = self.y_eng / b_half  # dimensionless lateral coordinate of engine attachment

        R_in_eng = 3 * ((eta_eng ** 2) / eta_cp) * (self.We / self.MTOM)  # engine relief factor

        W_f = self.MTOM - self.W_Zf  # fuel weight

        eta_f = (1 + 2 * self.TR_ft + 3 * self.TR_ft ** 2) / (4 * (1 + self.TR_ft + self.TR_ft ** 2)) * self.b_ft / self.b

        M_LF_y0 = 0.25 * eta_cp * b_st * n_ult * self.MTOM * (
                1 - eta_f / eta_cp * (1 - (self.W_Zf / self.MTOM)))  # upward moment due to lift and fuel
        R_in_f = 0.5 * (self.b_ft / self.b) * (1 + ((3 * self.TR ** 2) / (1 + 2 * self.TR))) * (1 - (self.W_Zf / self.MTOM))

        delta_nid = 10 ** -3  # smeared thickness increment to panel skins and webs
        S_box = self.S * 0.7  # total surface area of the wing box (box panels + spar webs)
        dW_nid_1 = 1.2 * delta_nid * rho_g * self.S  # weight increment due to non-taper, joints, fasteners
        dW_nid_2 = 0.18 * n_ult * WG * eta_cp * b_st * (rho_g / sigma)  # weight increment due to fail safety and damage tolerance

        omega_ic = 0.25

        R_ic = 1 + 2 * ((omega_ic * b_st) / self.S)  # correction factor for manholes and access hatches

        dW_nid_3 = 0.015 * (1 + 0.2 * self.N_eng) * self.W_pp  # weight penalty for engine installation

        if self.choice == "Yes":
            dW_nid_4 = 0.0015 * self.n_land * self.W_ML
        else:
            dW_nid_4 = 0

        dW_nid_5 = 0.0003 * n_ult * self.MTOM

        gamma = 1.4
        p = 58000  # pressure at flight altitude
        q_D = (gamma / 2 * p * self.M_D ** 2) / 1000  # dynamic pressure at design diving speed [kNm^-2]
        W_ref = 200  # calibration factor [N]
        b_ref = 50  # reference span [m]
        t_over_c_ref = 0.1  # thickness ratio at reference wing chord just inside the aileron
        q_ref = 30  # reference dynamic pressure [kNm^-2]
        lambda_05_rad = self.lambda_05 * (np.pi / 180)
        lambda_0_rad = self.lambda_0 * (np.pi / 180)

        dW_nid_6 = W_ref * (q_D / q_ref) * ((self.b * np.cos(lambda_0_rad) / b_ref) ** 3) * (t_over_c_ref) ** -2 * (1 - np.sin(lambda_05_rad)) * (1 - (self.M_D * np.cos(self.lambda_05)) ** 2) ** -0.5

        # Secondary Structures and Miscellaneous Items

        S_ref = 10  # reference area [m^2]
        omega_ref = 56  # reference specific weight [N/m^2]
        W_ref_1 = 10 ** 6 # reference weight [N]
        if self.LE_slat_choice == "Yes":
            k_fle = 1.3
        else:
            k_fle = 1

        omega_fle = 3.15 * k_fle * omega_ref * (q_D / q_ref) ** 0.25 * (
                    (MTOW * b_st) / (W_ref_1 * b_ref)) ** 0.145  # specific weight of fixed leading edge [N/m^2]
        omega_slat = 4.83 * omega_ref * (self.S_slat / S_ref) ** 0.183  # specific weight of LE high-lift devices [N/m^2]
        omega_fte = 2.6 * omega_ref * ((MTOW * b_st) / (W_ref_1 * b_ref)) ** 0.0544

        if self.TE_flap_choice == "Single or Double Slotted Fowler Flaps":
            omega_fte += 40
            k_slot = 1
        else:
            omega_fte += 100
            k_slot = 1.5

        k_sup = 1

        if self.ail_choice == "Unbalanced Ailerons":
            k_bal = 1

        elif self.ail_choice == "Balanced Ailerons":
            k_bal = 1.3

        else:
            k_bal = 1.54

        omega_tef = 1.7 * k_sup * k_slot * omega_ref * (
                    1 + (MTOW / W_ref_1) ** 0.35)  # specific weight of trailing edge flaps [N/m^2]

        omega_ail = 3 * omega_ref * k_bal * (self.S_ail / S_ref) ** 0.044  # aileron specific weight [N/m^2]
        omega_sp = 2.2 * omega_ref * (self.S_sp / S_ref) ** 0.032  # spoiler specific weight [N/m^2]

        W_tip = 150 * (MTOW / W_ref_1) * 0.67  # [tip structure weight]

        l_ref = 5
        omega_tip = 2.5 * omega_ref * (
                    (MTOW * self.l_tip) / (W_ref_1 * l_ref)) ** 0.145  # wingtip (winglet) specific weight [N/m^2]

        # Stress levels in aluminium alloys

        sigma_y_T3 = 450  # yield stress of 2024-T3 [MPa]
        sigma_y_7075_T6 = 520  # yield stress of 7075-T6 [MPa]
        sigma_y_7150_T6 = 550  # yield stress of 7150-T6 [MPa]
        sigma_y_T7 = 620  # yield stress of 7055-T7 [MPa]

        R_in = 1 - (R_in_f + R_in_wg + R_in_eng)
        SI = (0.5 * n_ult * MTOW * ((R_in * eta_cp * b_st) / (
                    eta_t * t_r * self.c_box * self.r))) / 10 ** 6  # structural index in compression [MPa]

        tau_mean = 0.5 * sigma

        # Refinements

        R_tip = 1 + (self.S_tip / self.S) * (((1 - eta_cp) / eta_cp) * np.cos(gamma_tip_rad) + (
                    self.b_tip / self.b))  # correction factor on cantilever ratio due to tip extension
        R_cs = (1 - (self.b_cs / self.b)) * (1 + (self.b_cs / self.b) * (2.45 * np.cos(lambda_ea) - 1))
        R_prime_cant = R_tip * ((self.b - self.b_cs) / (2 * t_r * np.cos(lambda_ea))) * (2 / 3 + ((t_r / self.c0) / (3 * (self.t_tk / self.c_tk))))  # modification for cantilever ratio

        # Advanced Materials

        if self.choice_fle == "Yes":
            R_fle = 0.8  # weight reduction due to CFC for fixed leading edge
        else:
            R_fle = 1
        if self.choice_tef == "Yes":
            R_tef = 0.85  # weight reduction due to CFC for trailing edge flaps
        else:
            R_tef = 1
        if self.choice_ail_cfc == "Yes":
            R_ail = 0.85
        else:
            R_ail = 1
        if self.choice_spoil == "Yes":
            R_sp = 0.9  # weight reduction due to CFC for spoilers
        else:
            R_sp = 1
        if self.choice_w_fus == "Yes":
            R_w_fus = 0.8  # weight reduction due to CFC for wing/fuselage fairing
        else:
            R_w_fus = 1

        omega_fle = omega_fle * R_fle
        omega_tef = omega_tef * R_tef
        omega_ail = omega_ail * R_ail
        omega_sp = omega_sp * R_sp

        # Application

        sigma_r_appl = 1 / (0.5 * ((R_ic / sigma_t) + (1.25 / sigma_c)))
        W_id_box = 0.36 * n_ult * R_in * WG * eta_cp * b_st * (rho_g / sigma_r_appl) * (1.05 * (R_cant / eta_t) + 3.67)
        W_id = W_id_box + W_rib


        W_id_box_prime = 0.36 * n_ult * R_tip * R_cs * WG * eta_cp * b_st * (rho_g / sigma_r_appl) * (
                1.05 * (R_prime_cant / eta_t) + 3.67)
        W_id_prime = (W_id_box_prime + W_rib)/self.g

        W_fle = self.S_fle * omega_fle
        W_fte = self.S_fte * omega_fte
        W_slat = self.S_slat * omega_slat
        W_sp = self.S_sp * omega_sp
        W_tip = self.S_tip * omega_tip
        W_ail = self.S_ail * omega_ail

        W_sec = (W_fle + W_fte + W_slat + W_sp + W_tip + W_ail)/self.g

        W_nid = (dW_nid_1 + dW_nid_2 + dW_nid_3 + dW_nid_4 + dW_nid_5 + dW_nid_6)/self.g

        W_wing = W_id_prime + W_nid + 1.1 * (W_sec)

        results_text += f"Ideal Weight (kg): {W_id_prime}\n"
        results_text += f"Non-Ideal Weight (kg): {W_nid}\n"
        results_text += f"Weight of Secondary Items (kg): {W_sec}\n"
        results_text += f"Total Wing Weight (kg): {W_wing}\n"

        return results_text
