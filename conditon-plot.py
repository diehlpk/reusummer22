import numpy as np
import matplotlib.pyplot as plt

pgf_with_latex = {"text.usetex": True, "font.size" : 12, "pgf.preamble" : [r'\usepackage{xfrac}'] }

# Mixed 

con_vhcm_n = np.genfromtxt('con_vhcm_n-damage.csv',delimiter=',')
con_mscm_n = np.genfromtxt('con_mscm_n-damage.csv',delimiter=',')
con_mdcm_n = np.genfromtxt('con_mdcm_n-damage.csv',delimiter=',')

x = np.linspace(0,1,len(con_vhcm_n))



plt.plot(x,con_mdcm_n,label="MDCM",c="black",marker="s")
plt.plot(x,con_mscm_n,label="MSCM",c="black",marker="o")
plt.plot(x,con_vhcm_n,label="VHCM",c="black",marker="x")
plt.grid()
plt.xlabel("$\delta$")
plt.ylabel("Condition number")
plt.yscale('log')
plt.legend()
plt.title("Mixed boundary conditions \n Quadratic interpolation")
ax = plt.gca()
ax.set_xticks(x)
ax.set_xticklabels(["1/8","1/16","1/32","1/64"])
plt.savefig("Condition-n-damage.pdf",bbox_inches='tight')
