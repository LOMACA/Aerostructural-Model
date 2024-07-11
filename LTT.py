import numpy as np

class LiftingLineTheory:
    def __init__(self, v_inf, c_root, c_tip, b, n, cl0, cd0, a0, alpha0, rho):
        self.v_inf = v_inf
        self.c_root = c_root
        self.c_tip = c_tip
        self.b = b
        self.n = n
        self.cl0 = cl0
        self.cd0 = cd0
        self.a0 = a0
        self.alpha0 = alpha0 * np.pi / 180
        self.rho = rho


    def compute(self, alpha_values):

        results = []

        b_half = self.b / 2
        TR = self.c_tip / self.c_root
        Aw = b_half * (self.c_root + self.c_tip)
        S = 2 * Aw
        AR = self.b ** 2 / (2 * Aw)

        if self.n == 4:
            phi = np.array([22.5, 45, 67.5, 90])
        else:
            phi = np.linspace(0.1e-8, 90, self.n)

        phi_r = phi * np.pi / 180

        y = (self.b / 2) * np.cos(phi_r)

        c = self.c_root * (1 + (TR - 1) * np.cos(phi_r))
        geo = (c * self.a0) / (2 * AR * self.c_root * (1 + TR))

        for alpha in alpha_values:
            alpha_rad = alpha * np.pi / 180

            Left = geo * (alpha_rad - self.alpha0) * np.sin(phi_r)
            B = np.zeros((self.n, self.n))

            for i in range(self.n):
                for j in range(self.n):
                    B[i, j] = np.sin((2 * (j + 1) - 1) * phi_r[i]) * (geo[i] * (2 * (j + 1) - 1) + np.sin(phi_r[i]))

            A = np.linalg.solve(B, Left)

            CL = np.pi * AR * A[0]
            delta = np.sum([(2 * (i + 1) - 1) * A[i] ** 2 for i in range(1, self.n)])
            CD_ind = (CL ** 2 / (np.pi * AR)) * (1 + delta)

            CD = self.cd0 + CD_ind

            L = CL * 0.5 * self.rho * self.v_inf ** 2 * Aw * 2  # Lift [N]
            x = np.linspace(-self.b / 2, self.b / 2, self.n)
            D_ind = CD_ind * 0.5 * self.rho * self.v_inf ** 2 * Aw  # Induced Drag [N]
            D = CD * 0.5 * self.rho * self.v_inf ** 2 * Aw  # Total Drag [N] (neglecting wave drag and viscous drag)
            gamma_temp = np.zeros((self.n))
            for i in range(0, self.n):
                for j in range(0, self.n):
                    gamma_temp[i] = gamma_temp[i] + A[j] * np.sin((2 * (j + 1) - 1) * phi_r[i])

            gamma = 2 * self.b * self.v_inf * gamma_temp  # Circulation distribution [m^2/s]
            tau = 0.05
            aw = self.a0 / (1 + (self.a0 / (np.pi * AR)) * (1 + tau))  # lift slope of finite wing [rad]
            L_dist = self.rho * self.v_inf * gamma  # lift distribution [N]
            M_sum = np.zeros_like(L_dist)
            y_new = y[::-1]
            L_dist_new = L_dist[::-1]  # re-ordered lift distribution [N]
            for i in range(len(y)):
                if i == 0:
                    continue  # Skip the first element (tip of the wing)
                # Sum the contributions from the current point to the wing root
                M_sum[i] = np.sum(L_dist_new[:i] * np.diff(y_new[:i + 1]))
            M_dist = M_sum  # bending moment distribution [Nm]

            temp = np.zeros(self.n)
            for i in range(0, self.n):
                for j in range(0, self.n):
                    temp[i] = temp[i] + ((2 * (j + 1) - 1)) * A[j] * np.sin((2 * (j + 1) - 1) * phi_r[i]) / np.sin(
                        phi_r[i])

            dw = -self.v_inf * temp  # downwash [m/s]
            alpha_i = dw / self.v_inf  # induced angle of attack [rad]

            D_dist = L_dist * np.abs(alpha_i)  # Induced Drag [N]

            results.append({
                'alpha': alpha,
                'CL': CL,
                'CD': CD,
                'Lift': CL * 0.5 * self.rho * self.v_inf ** 2 * Aw * 2,
                'Drag': CD * 0.5 * self.rho * self.v_inf ** 2 * Aw * 2,
                'Gamma': gamma,
                'M_dist': M_dist,
                'D_dist': D_dist,
                'L_dist': L_dist,
                'dw': dw,
                'y': y,
                'aw': aw,
                'AR': AR,
                'TR': TR,
                'S': S,
                'v_inf': self.v_inf,
                'c_root': self.c_root,
                'c_tip': self.c_tip,
                'b': self.b,
                'n': self.n
            })

        return results

