#!/usr/bin/env python3
# Coupling using the displacement for MDCM
# @author patrickdiehl@lsu.edu
# @author serge.prudhomme@polymtl.ca
# @author a30110572@mcla.edu
# @date 02/05/2021
import numpy as np
import sys 
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.ticker import FormatStrFormatter

# pgf_with_latex = {"text.usetex": True, "font.size" : 12, "pgf.preamble" : [r'\usepackage{xfrac}'] }


example = sys.argv[1]

g = -1
has_condition = True
con = []

#############################################################################
# Solve the system
#############################################################################

def solve(M,f):
    
    return np.linalg.solve(M,f)

#############################################################################
# Loading
#############################################################################

def f(x):
    
    global g 

    if example == "Cubic":
        g = 27
        return -6*x
    elif example == "Quartic":
        g = 108
        return -12 *x*x
    elif example == "Quadratic":
        g = 6
        return -2
    elif example == "Linear":
        g = 1
        return 0.1
    elif example == "Linear-cubic":
        g = 31./4.
        if x < 1.5:
            return 0.1
        else:
            return 9-6*x
    else:
        print("Error: Either provide Linear, Quadratic, Quartic, or Cubic")
        sys.exit()

def forceFull(n,h):
    
    force = np.zeros(n)
   
    for i in range(1,n-1):
        force[i] = f(i * h)
    
    force[n-1] = g
    
    return force

def forceCoupling(n,x):
    
    force = np.zeros(3*n+4)
   
    for i in range(1,3*n+3):
        force[i] = f(x[i])
    
    force[3*n+3] = g
    
    return force

def forceCouplingFD(n,x):
    
    force = np.zeros(3*n+3)
   
    for i in range(1,3*n+2):
        force[i] = f(x[i])
    
    force[3*n+2] = g
    
    return force

#############################################################################
# Exact solution 
#############################################################################

def exactSolution(x):
    
    if example == "Cubic":
        return x * x * x
    elif example == "Quartic":
        return x * x * x * x
    elif example == "Quadratic":
        return x * x
    elif example == "Linear":
        return x
    elif example == "Linear-cubic":
        return np.where(x < 1.5, x, x + (x-1.5) * (x-1.5) * (x-1.5) )
    else:
        print("Error: Either provide Linear, Quadratic, Quartic, or Cubic")
        sys.exit()

#############################################################################
# Assemble the stiffness matrix for the finite difference model (FD)
#############################################################################

def FDM(n,h):

    M = np.zeros([n,n])

    M[0][0] = 1

    for i in range(1,n-1):
        M[i][i-1] = -2
        M[i][i] = 4
        M[i][i+1] = -2

    M[n-1][n-1] = 11*h / 3
    M[n-1][n-2] = -18*h / 3
    M[n-1][n-3] = 9 * h / 3
    M[n-1][n-4] = -2 * h / 3


    M *= 1./(2.*h*h)

    return M

#############################################################################
# Assemble the stiffness matrix for the coupling of FDM - FDM - FDM
#############################################################################

def CouplingFDFD(n,h):

    M = np.zeros([3*n,3*n])

    fFD = 1./(2.*h*h)

    M[0][0] = 1

    for i in range(1,n-1):
        M[i][i-1] = -2 * fFD
        M[i][i] = 4 * fFD
        M[i][i+1] = -2 * fFD

    M[n-1][n-1] = -1 
    M[n-1][n] = 1 

    M[n][n-1] = 3*h * fFD
    M[n][n-2] = -4*h * fFD
    M[n][n-3] = 1*h * fFD

    M[n][n] = 3*h * fFD
    M[n][n+1] = -4*h * fFD
    M[n][n+2] = 1*h * fFD

    for i in range(n+1,2*n-1):
        M[i][i-1] = -1 * (((fPD(x[i],h)) + (fPD(x[i-1],h)))/2)
        M[i][i] = 1 * (2 * (fPD(x[i],h)))
        M[i][i+1] = -1 * (((fPD(x[i],h)) + (fPD(x[i+1],h)))/2)

    M[2*n-1][2*n-1] = -1 
    M[2*n-1][2*n] = 1

    M[2*n][2*n-1] = 3*h * fFD
    M[2*n][2*n-2] = -4*h * fFD
    M[2*n][2*n-3] = h * fFD

    M[2*n][2*n] = 3*h * fFD
    M[2*n][2*n+1] = -4*h * fFD
    M[2*n][2*n+2] = h * fFD

    for i in range(2*n+1,3*n-1):
        M[i][i-1] = -2 * fFD
        M[i][i] = 4 * fFD
        M[i][i+1] = -2 * fFD

    M[3*n-1][3*n-1] = 3*h * fFD
    M[3*n-1][3*n-2] = -4*h * fFD
    M[3*n-1][3*n-3] = h * fFD

    #M *= 1./(2.*h*h)
    
    return M

#############################################################################
# Assemble the stiffness matrix for the coupling of FDM - Displacement - FDM 
#############################################################################

c = 0.5

def fPD(x,h):
    E = 1
    if x >= 1.25 and x <= 1.5:
        E = 1+4*(1-c)*(1.25-x)
    elif x >= 1.5 and x <= 1.75:
        E = 1+4*(1-c)*(x-1.75)
    return E/(h*h)
    

def Coupling(n,h,x):

    M = np.zeros([3*n+4,3*n+4])

    fFD =  1./(2.*h*h)
    # fPD =  1./(8.*h*h)

    # Boundary

    M[0][0] = 1

    # FD 

    for i in range(1,n-1):
        M[i][i-1] = -2 * fFD
        M[i][i] = 4 * fFD
        M[i][i+1] = -2 * fFD

    # Overlapp

    M[n-1][n-1] = -1
    M[n-1][n+2] = 1

    M[n][n] = -1
    M[n][n-3] = 1

    M[n+1][n+1] = -1
    M[n+1][n-2] = 1

    # PD

    for i in range(n+2,2*n+2):
        M[i][i-2] = 1.  * (fPD(x[i-2],h)/4)
        M[i][i-1] = 1. * (fPD(x[i-1],h))
        M[i][i] = -1 * (fPD(x[i-2],h)/4 + (fPD(x[i-1],h)) + (fPD(x[i+1],h)) + fPD(x[i+2],h)/4)
        M[i][i+1] =  1. * (fPD(x[i+1],h))
        M[i][i+2] = 1. * (fPD(x[i+2],h)/4)

    # Overlap

    M[2*n+2][2*n+2] = -1
    M[2*n+2][2*n+5] = 1

    M[2*n+3][2*n+3] = -1
    M[2*n+3][2*n+6] = 1

    M[2*n+4][2*n+4] = -1
    M[2*n+4][2*n+1] = 1

    # FD

    for i in range(2*n+5,3*n+3):
        M[i][i-1] = -2 * fFD
        M[i][i] = 4 * fFD
        M[i][i+1] = -2 * fFD

    # Boundary
 
    M[3*n+3][3*n+3] = 11 *  h * fFD / 3
    M[3*n+3][3*n+2] =  -18 * h * fFD  / 3
    M[3*n+3][3*n+1] = 9 * h * fFD / 3
    M[3*n+3][3*n] = -2 * h * fFD / 3

    if has_condition:
        con.append(np.linalg.cond(M))

    return M


markers = ['s','o','x','.']

for i in range(4,8):
    n = np.power(2,i)
    h = 1./n
    nodes = n + 1
    nodesFull = 3 * n + 1

    print(nodes,h)
    x1 = np.linspace(0,1,nodes)
    x2 = np.linspace(1-2*h,2+2*h,nodes+4)
    x3 = np.linspace(2,3.,nodes)
    x = np.array(np.concatenate((x1,x2,x3)))
    x1FD = np.linspace(0,1,nodes)
    x2FD = np.linspace(1,2,nodes)
    x3FD = np.linspace(2,3.,nodes)
    xFD = np.array(np.concatenate((x1FD,x2FD,x3FD)))

    xFull = np.linspace(0,3.,nodesFull)
    

  
    forceCoupled = forceCoupling(nodes,x)

    forceCoupled[nodes-1] = 0
    forceCoupled[nodes] = 0
    forceCoupled[nodes+1] = 0

    forceCoupled[2*nodes+2] = 0
    forceCoupled[2*nodes+3] = 0
    forceCoupled[2*nodes+4] = 0

    forceCoupledFD = forceCouplingFD(n,xFD)

    forceCoupledFD[n] = 0
    forceCoupledFD[2*n+1] = 0

    uFDMVHM = solve(Coupling(nodes,h,x),forceCoupled)
    uFD = solve(CouplingFDFD(nodes,h),forceCoupledFD)

    uSlice = np.array(np.concatenate((uFDMVHM[0:nodes],uFDMVHM[nodes+3:2*nodes+2],uFDMVHM[2*nodes+5:len(x)])))
    uSliceFD = np.array(np.concatenate((uFD[0:nodes],uFD[nodes+1:2*nodes],uFD[2*nodes+1:len(x)])))
    

    plt.axvline(x=1,c="#536872")
    plt.axvline(x=2,c="#536872")
    
    if example == "Quartic" or example == "Linear-cubic" or example =="Linear" or example == "Cubic" or example == "Quadratic":
        
        plt.plot(xFull,uSlice,label=r"$\delta$=1/"+str(int(n/2))+"",c="black",marker=markers[i-4],markevery=n)
        plt.plot(xFull,uSliceFD,label=r"$\delta$=1/"+str(int(n/2))+"",c="red",marker=markers[i-8],markevery=n)
        plt.ylabel("Error in displacement w.r.t. FDM")
        plt.ylabel("Error in displacement w.r.t. FDM")

    elif i == 4:

        plt.plot(xFull,uFD,label="FDM",c="black")
        plt.plot(xFull,uSlice,label=r"$\delta$=1/"+str(int(n/2))+"",c="black",marker=markers[i-4],markevery=n)
        plt.ylabel("Displacement")
        np.savetxt("coupling-"+example.lower()+"-approach-1.csv",uSlice)   

# plt.gca().yaxis.set_major_formatter(FormatStrFormatter('%0.5f')) 
plt.title("Example with "+example.lower()+" solution for MDCM with $m=2$")
plt.legend()
plt.grid()
plt.xlabel("$x$")

plt.savefig("coupling-"+example.lower()+"-approach-1.pdf",bbox_inches='tight')

if has_condition :
    file = "con_mdcm_d-cubic-matching.csv" # + str(c) + ".csv"
    # print(file)
    np.savetxt(file, con, delimiter=",")

