# PEM Electrolyzer model v0.1.0
from math import log10, exp
from matplotlib.pyplot import plot, show


# Model parameters to be included in config.ini
T = 25.0                                            # (°C) Electrolyzer operating temperature
Ir = 80.0                                           # (A) rated current
A = 0.25                                            # (m^2) area of electrode
Nc = 12.0                                           # Number of cells connected in series
F = 96485.34                                        # (C/mol) Faraday's constant
ne = 2.0
DG_25 = 237000.0
DG_80 = 228480.0
DG = DG_25-(T-25)/55*(DG_25-DG_80)
V_rev = DG/ne/F
DH = 286000.0
Vtn = DH/ne/F                                       # thermo-neutral voltage
r1 = 7.331e-5                                       # (ohm m^2) ri parameter for ohmic resistance of electrolyte
r2 = -1.107e-7                                      # (ohm m2 °C^-1)
r3 = 0
s1 = 1.586e-1                                       # (V) si and ti parameters for over-voltage on electrodes
s2 = 1.378e-3                                       # (V°C^-1)
s3 = -1.606e-5                                      # V °C^-2)
t1 = 1.599e-2
t2 = -1.302
t3 = 4.213e2
I_initial = 0
I_final = 870
I_step = 1
nn = (I_final-I_initial)/I_step
I = [i for i in range(int(nn) + 1)]

# Compute V and Id
V = [V_rev + (r1+r2*T)*I[i]/A+(s1+s2*T+s3*T**2)*log10((t1+t2/T+t3/T**2)*I[i]/A+1) for i in range(int(nn) + 1)]
Id = [i/2500 for i in I]

# Compute rated voltage (Vr)
Vr = V_rev + (r1+r2*T)*Ir/A+(s1+s2*T+s3*T**2)*log10((t1+t2/T+t3/T**2)*Ir/A+1)
P = Nc*Vr*Ir
plot(Id, V)
show()

# Faraday Efficiency
a1 = 0.995                                          # 99.5 %
a2 = -9.5788                                        # (m ^ 2 * A ^ -1)
a3 = -0.0555                                        # (m ^ 2 * A ^ -1 *°C)
a4 = 0
a5 = 1502.7083                                      # (m ^ 4 * A ^ -1)
a6 = -70.8005                                       # (m ^ 4 * A ^ -1 *°C-1)
a7 = 0
nf = a1*exp((a2+a3*T+a4*T**2)/(Ir/A)+(a5+a6*T+a7*T**2)/(Ir/A)**2)

# Energy (or voltaje) efficiency of a cell
ne = Vtn/Vr

# Flow of H2 produced
Qh2 = 80.69*Nc*Ir*nf/2/F                            # (Nm^3/h)
