function v_eq = ss_equilibrium(g_m, g_n, I)

v_eq = (I + 40*g_n - 80*g_m)./(g_n + g_m);