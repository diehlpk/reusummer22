import numpy as np
import matplotlib.pyplot as plt

pgf_with_latex = {"text.usetex": True, "font.size" : 12, "pgf.preamble" : [r'\usepackage{xfrac}'] }

# Mixed 

c = [0.1,0.25,0.75,0.9]

x = np.linspace(0,1,len(c))
markers = ['s','o','x','.']

i = 0
for value in c:
    con = np.genfromtxt("con_neumann_"+str(value)+".csv",delimiter=',')
    plt.plot(x,con,label=str(value),c="black",marker=markers[i])
    i += 1

plt.grid()
plt.xlabel("$\delta$")
plt.ylabel("Condition number")
plt.yscale('log')
plt.legend()
plt.title("Mixed boundary conditions \n Damage")
ax = plt.gca()
ax.set_xticks(x)
ax.set_xticklabels(["1/8","1/16","1/32","1/64"])
plt.savefig("Condition-n-damage.pdf",bbox_inches='tight')

plt.clf()

i = 0
for value in c:
    con = np.genfromtxt("con_d_"+str(value)+".csv",delimiter=',')
    print(con)
    plt.plot(x,con,label=str(value),c="black",marker=markers[i])
    i += 1

plt.grid()
plt.xlabel("$\delta$")
plt.ylabel("Condition number")
plt.yscale('log')
plt.legend()
plt.title("Homogeneous boundary conditions \n Damage")
ax = plt.gca()
ax.set_xticks(x)
ax.set_xticklabels(["1/8","1/16","1/32","1/64"])
plt.savefig("Condition-d-damage.pdf",bbox_inches='tight')
