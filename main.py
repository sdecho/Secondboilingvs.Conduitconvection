import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# A conductive cooling model by Annen et al.(2006)

density_basalt = 2480  # 2672 # kg/m3
density_uc = 3050 # kg/m3
cp_uc = 1370 # J/(kg*K)
cp_basalt = 1480 # J/(kg*K)
L_basalt = 4*(10**5)
L_uc = 2.7*(10**5)
k0_basalt = 2.6 #J/(s.K)
k0_uc = 3 # J/(s.k)


Tb_i = 1150 + 273
Ts = 930+273
Tl = 1145 +273


# create two vectors for time and depth from the top to the center of magma chamber
tmax = 1000*365*24*3600
n_t = 1000
t = np.linspace(0,tmax, num=n_t)
dt = t[1]-t[0]

n_d = 1000
dmax = 10000 #m
d = np.linspace(0, dmax, num=n_d)
dd = d[1]-d[0]

# initialize a matrix for temperature with n_t rows and n_d columns
initial = 400
T = np.zeros((n_t,n_d))
fluid = np.zeros((n_t,n_d))
T[:,0] = 273 + initial

D_up = 4000
D_low = 7000
gradient = 40

# define T profile at t = 0 along depth
for i in range (0, n_d):
    if d[i] < D_up:
        temp = initial+gradient*d[i]/1000
        T[0, i] = temp +273 # in K
    elif d[i] <= D_low:
        T[0, i] = Tb_i
    else:
        T[0, i] = initial+273 + gradient*d[i]/1000
T[:, n_d-1] = 10*gradient+273+initial
index_c = np.where(d > D_up+150)
center = index_c[0][0]
index_up = np.where(d > D_up)
up = index_up[0][0]
Tc = [T[0, center]-273]
Tup = [T[0, up]-273]

H2O_solubility = 2.5
S_fluid = 0.15
S_total = 0
S_flux = [0]
S_flux_2 = [0]
print(dd)
# solve the equation
for p in range(0, n_t-1):
    S_pertime = 0
    for i in range(1, n_d-1):
        if d[i] < D_up and d[i] > D_low:
            A = 0
            k = k0_uc * (1 + 1.5 * (10 ** -3) * d[i] / 1000) / (1 + T[p, i] * (10 ** -4))
            L = L_uc
            density = density_uc
            cp = cp_uc
        else:
            A = 3.25 * (10**-3)
            k = k0_basalt * (1 + 1.5 * (10 ** -3) * d[i] / 1000) / (1 + T[p, i] * (10 ** -4))
            if T[p, i] < Tl and T[p, i] > Ts:
                L = L_basalt
            else:
                L = 0
            density = density_basalt
            cp = cp_basalt
        T[p+1, i] = T[p, i] + k * dt/((density*cp+A*L)*(dd**2))*(T[p, i+1]-2*T[p,i]+T[p, i-1])
        if T[p+1, i] < Tl and T[p+1, i] > Ts:
            if (-0.0049 * T[p+1, i] + 6.8904) > 0:
                crystal = -0.0049 * T[p+1, i] + 6.8904
            else:
                crystal = 0
            fluid[p+1, i] = H2O_solubility*crystal/100
            S_pertime = S_pertime + (fluid[p+1, i]-fluid[p, i])*dd*(3000*3000)*density_basalt * S_fluid
    S_flux_2.append(S_pertime/(dt/(24*3600)))
    S_total = S_total + S_pertime
    S_avflux = (S_total/1000)/(t[p+1]/(24*3600))
    S_flux.append(S_avflux)

    # S_pertime =S_pertime + S_total
    # S_avflux = S_pertime/(t[p+1]/(3600*24))/1000
    Tc.append(T[p+1, center]-273)
    Tup.append(T[p+1, up]-273)

        # if d[i] == dmax:
        #     T[p + 1, i] = T[p, i] + k_uc * dt / ((density_uc * cp_uc + A * L_uc) * (dd ** 2)) * (T[p, i - 1] - T[p, i])
S_averageflux = sum(S_flux)/len(S_flux)


print(f"The average S flux is {S_averageflux} T/day")


###########
# This calculation is derived from Stevenson and Blake, 1998
density_c = 2410 # TRpreRHcr1_mi1, 1160 C, 2.5 wt.%
density_d = 2507 # 2953
viscosity_c = np.exp(2.3786313019514123)*((1-0.05/0.66)**-2)  # Pa/s 5vol% crystals
# viscosity_d = np.exp(5.483709785192866)*((1-0.34/0.66)**-2) # Pa/s 34 vol% crystals
viscosity_d = np.exp(6.07)*((1-0.34/0.66)**-2) # Pa/s 34 vol% crystals
# viscosity_d = np.exp(6)
print(viscosity_d, viscosity_c, viscosity_d/viscosity_c)


Ps = 0.064 #0.02+0.04*np.log10(viscosity_d/viscosity_c)
R_start = 0.6
R = np.linspace(0.01,5, num=100)

Q_m = 3.14 *(R_start)**2*Ps*(9.8*(density_d-density_c)*R**4)/viscosity_d
v = Q_m/(3.14*(R_start*R)**2)
print(v)
dc_S = 0.15
Q_S = Q_m * density_c * dc_S*0.01 * 3600*24/1000
v_mb_1 = 2.2*250*1000/(density_c*3.14*(R_start*R)**2*dc_S*0.01*24*3600)
v_mb_2 = 2.2*375*1000/(density_c*3.14*(R_start*R)**2*dc_S*0.01*24*3600)
print(2.2*250*1000/(density_c*3.14*(R_start*5)**2*dc_S*0.01*24*3600))

for r in R:
    q = 3.14 * (R_start) ** 2 * Ps * (9.8 * (density_d - density_c) * r ** 4) / viscosity_d
    velo = q/(3.14*(R_start*r)**2)
    v_250 = 2.2*250 * 1000 / (density_c * 3.14 * (R_start * r) ** 2 * dc_S * 0.01 * 24 * 3600)
    v_375 = 2.2*375 * 1000 / (density_c * 3.14 * (R_start * r) ** 2 * dc_S * 0.01 * 24 * 3600)

    if abs(velo- v_250)<0.0001:
        R_250 = r

    elif abs(velo- v_375)<0.0001:
        R_375 = r
        print(v_375, R_375)

R_250 = 2.6200
v_250= 0.2269

R_375= 2.8992
v_375 = 0.2779

#############################################
# Figure 9
fig3 = plt.figure()
ax1 = fig3.add_subplot(2,2,1)
ax1.plot(T[0, :] - 273 * np.ones(n_d), d, "black", linestyle=":", linewidth=4)
ax1.plot(T[n_t-1,:] - 273 * np.ones(n_d), d, "black", linestyle="-")
ax1.set_xlabel("Temperature ($^\circ$C)")
ax1.set_ylabel("depth (km)")
ax1.set_ylim([10000, 0])
ax1.set_xlim([0, 1200])
ax1.tick_params(axis="x", colors="black")
ax1.legend(["$Temperature_{time=0}$", "$Temperature_{time=1000y}$"], loc="upper right", frameon=False, title="temperature profile")
ax1.annotate("(a)", xy=(0.02, 0.90), xycoords="axes fraction", fontsize=14, weight='regular', color="black")
ax2 = fig3.add_subplot(221, frame_on=False)
ax2.plot(fluid[n_t-1, :]*100, d, "grey", linestyle="-", linewidth=3)
ax2.xaxis.tick_top()
ax2.xaxis.set_label_position('top')
ax2.tick_params(axis='x', colors="grey")
ax2.set_ylim([10000,0])
ax2.set_xlabel("fluid (wt.%)")
ax2.xaxis.label.set_color('grey')
ax2.set_xlim([2, 0])
ax2.legend(["$fluid_{time=1000y}$"], loc="lower right", frameon=False, title="fluid profile", labelcolor= "grey")
ax2.get_legend().get_title().set_color('grey')

ax3 = fig3.add_subplot(2,2,3)
ax3.plot(t[1:]/(3600*24*365), S_flux[1:],"blue")
ax3.set_xlim([0, 1000])
ax3.set_ylim([0, 1001])
ax3.set_yticks(np.arange(0, 1001, 100))
ax3.set_xlabel("time (year)")
ax3.set_ylabel("S flux (t/day)")
ax3.annotate(f"The average S flux is {round(S_averageflux,1)} t/day,\ndt = {round(dt/(365*24*3600),1)} year",
             xy=(0.05, 0.85), xycoords="axes fraction", fontsize=14, weight='regular', color="blue")
ax3.annotate("(b)", xy=(0.02, 0.05), xycoords="axes fraction", fontsize=14, weight='regular', color="black")

ax4 = fig3.add_subplot(2,2,4)
ax4.plot(R, Q_S, "blue", linestyle="--" )
ax4.legend(["S&B model"])
ax4.set_ylim([0,4000])
ax4.set_xlim([0,5])
ax4.set_xlabel("conduit radius (m)")
ax4.set_ylabel("S flux (t/day)")
ax4.annotate("(d)", xy=(0.90, 0.05), xycoords="axes fraction", fontsize=14, weight='regular', color="black")

ax5 = fig3.add_subplot(2,2,2)
ax5.plot(R, v, "blue", linestyle="--")
ax5.plot(R, v_mb_1, "green")
ax5.plot(R, v_mb_2, "red")
ax5.plot(R_250, v_250, "*", markerfacecolor = "green",markersize = 15)
ax5.plot(R_375, v_375, "*", markerfacecolor = "red", markersize = 15)
ax5.legend(["S&B model", "S flux = 550 t/day", "S flux = 825 t/day"],loc='upper right')
ax5.set_xlim([0,5])
ax5.set_ylim([0,1])
ax5.set_xlabel("conduit radius (m)")
ax5.set_ylabel("ascent velocity (m/s)")
ax5.annotate("(c)", xy=(0.05, 0.90), xycoords="axes fraction", fontsize=14, weight='regular', color="black")


plt.show()
