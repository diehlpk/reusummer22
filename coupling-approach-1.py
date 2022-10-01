#!/usr/bin/env python3
# Coupling using the displacement for MDCM
# @author patrickdiehl@lsu.edu
# @author serge.prudhomme@polymtl.ca
# @author a30110572@mcla.edu
# @date 02/05/2021
from configparser import Interpolation
import numpy as np
import sys 
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.ticker import FormatStrFormatter
from scipy import interpolate

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
        g = 3/27*3*3
        return -6/27*x
    elif example == "Quartic":
        g = 4/81 * 3 * 3 * 3
        return -12/81 * x*x
    elif example == "Quadratic":
        g = 6/9
        return -2/9
    elif example == "Linear":
        g = 1/3
        return 0
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
        force[i] =  f(x[i])
    
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
        return 1/27 * x * x * x
    elif example == "Quartic":
        return x * x * x * x / 81
    elif example == "Quadratic":
        return 1/9 * x * x
    elif example == "Linear":
        return x/3
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

def CouplingFDFD(n,h,x):

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
        left =  0.5 * (fPD(x[i-1],h)  + fPD(x[i],h)  ) 
        right = 0.5 * (fPD(x[i],h)  + fPD(x[i+1],h) ) 
        M[i][i-1] =  - left
        M[i][i] =  (left + right)
        M[i][i+1] = - right

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

    M[3*n-1][3*n-1] = 11*h / 3 * fFD
    M[3*n-1][3*n-2] = -18*h / 3  * fFD
    M[3*n-1][3*n-3] = 9 * h / 3  * fFD
    M[3*n-1][3*n-4] =  -2 * h / 3  * fFD

    #M *= 1./(2.*h*h)

    #np.savetxt("foo.csv", M, delimiter=",")
    
    return M

#############################################################################
# Assemble the stiffness matrix for the coupling of FDM - Displacement - FDM 
#############################################################################

c = 0.1
x = [1.25,(1.5+1.25)/2,1.5,(1.5+1.75)/2,1.75]
y = [1,(1+c)/2,c,(1+c)/2,1]
tck = interpolate.splrep(x, y, s=0)

def fPD(x,h):
    E = 1
    if x >= 1.25 and x <= 1.75:
        E = interpolate.splev(x, tck, der=0)
    return E/(h*h)



def Coupling(n,h,x):

    M = np.zeros([3*n+4,3*n+4])

    fFD =  1./(2.*h*h)
    f =  1./(8.*h*h)

    # Boundary

    M[0][0] = 1

    # FD 

    for i in range(1,n-1):
        #print(x[i])
        M[i][i-1] = -2 * fFD
        M[i][i] = 4 * fFD
        M[i][i+1] = -2 * fFD

    # Overlapp
    #print("----")
    M[n-1][n-1] = -1
    M[n-1][n+2] = 1
    #print(x[n-1])
    #print(x[n+2]) 

    M[n][n] = -1
    M[n][n-3] = 1
    #print(x[n])
    #print(x[n-3])

    M[n+1][n+1] = -1
    M[n+1][n-2] = 1
    #print(x[n+1])
    #print(x[n-2])
    # PD
    #print("----")
    for i in range(n+2,2*n+2):
        if x[i] < 1.25 :
            M[i][i-2] = -1. * f
            M[i][i-1] = -4. * f
            M[i][i] = 10. * f
            M[i][i+1] = -4. * f
            M[i][i+2] = -1. * f
        elif x[i] > 1.75:
            M[i][i-2] = -1. * f
            M[i][i-1] = -4. * f
            M[i][i] = 10. * f
            M[i][i+1] = -4. * f
            M[i][i+2] = -1. * f
        #print(x[i])
        else:
            M[i][i-2] = -1.  * (0.5*(fPD(x[i-2],h)+fPD(x[i],h))/8) 
            M[i][i-1] = -1. * (0.5*(fPD(x[i-1],h)+fPD(x[i],h))/2) 
            M[i][i] = 1 * (0.5*(fPD(x[i-2],h)+fPD(x[i],h))/8+ 0.5*(fPD(x[i-1],h)+fPD(x[i],h))/2+ 0.5*(fPD(x[i+1],h)+fPD(x[i],h))/2  + 0.5*(fPD(x[i+2],h)+fPD(x[i],h))/8) 
            M[i][i+1] =  -1. * (0.5*(fPD(x[i+1],h)+fPD(x[i],h))/2) 
            M[i][i+2] = -1. * (0.5*(fPD(x[i+2],h)+fPD(x[i],h))/8) 

    # Overlap

    M[2*n+2][2*n+2] = -1
    M[2*n+2][2*n+5] = 1
    #print("----")
    #print(x[2*n+2])
    #print(x[2*n+5])

    M[2*n+3][2*n+3] = -1
    M[2*n+3][2*n+6] = 1
    #print(x[2*n+3])
    #print(x[2*n+6])

    M[2*n+4][2*n+4] = -1
    M[2*n+4][2*n+1] = 1
    #print(x[2*n+4])
    #print(x[2*n+1])
    #print("----")

    # FD

    for i in range(2*n+5,3*n+3):
        #print(x[i])
        M[i][i-1] = -2 * fFD
        M[i][i] = 4 * fFD
        M[i][i+1] = -2 * fFD

    # Boundary
    #print(x[3*n+3])
    M[3*n+3][3*n+3] = 11 *  h * fFD /3
    M[3*n+3][3*n+2] =  -18 * h * fFD /3  
    M[3*n+3][3*n+1] = 9 * h * fFD /3
    M[3*n+3][3*n] = -2 * h * fFD / 3

    if has_condition:
        con.append(np.linalg.cond(M))

    return M

def VHM(n,h,x):

    f = 1./(8.*h*h) 
    MVHM = np.zeros([n,n])

    MVHM[0][0] = 1. * f
    MVHM[1][0] = -8. * f
    MVHM[1][1] = 16. * f
    MVHM[1][2] = -8. * f

    for i in range(2,n-2):

        if x[i] < 1.25 :
            MVHM[i][i-2] = -1. * f
            MVHM[i][i-1] = -4. * f
            MVHM[i][i] = 10. * f
            MVHM[i][i+1] = -4. * f
            MVHM[i][i+2] = -1. * f
        elif x[i] > 1.75:
            MVHM[i][i-2] = -1. * f
            MVHM[i][i-1] = -4. * f
            MVHM[i][i] = 10. * f
            MVHM[i][i+1] = -4. * f
            MVHM[i][i+2] = -1. * f
        else:
            MVHM[i][i-2] = -1.  * (0.5*(fPD(x[i-2],h)+fPD(x[i],h))/8) 
            MVHM[i][i-1] = -1. * (0.5*(fPD(x[i-1],h)+fPD(x[i],h))/2) 
            MVHM[i][i] = 1 * (0.5*(fPD(x[i-2],h)+fPD(x[i],h))/8+ 0.5*(fPD(x[i-1],h)+fPD(x[i],h))/2+ 0.5*(fPD(x[i+1],h)+fPD(x[i],h))/2  + 0.5*(fPD(x[i+2],h)+fPD(x[i],h))/8) 
            MVHM[i][i+1] =  -1. * (0.5*(fPD(x[i+1],h)+fPD(x[i],h))/2) 
            MVHM[i][i+2] = -1. * (0.5*(fPD(x[i+2],h)+fPD(x[i],h))/8) 


    MVHM[n-2][n-1] = -8. * f
    MVHM[n-2][n-2] = 16. * f
    MVHM[n-2][n-3] = -8. * f

    MVHM[n-1][n-1] = 12.*h * f
    MVHM[n-1][n-2] = -16.*h * f
    MVHM[n-1][n-3] = 4.*h * f

    #MVHM *= 1./(8.*h*h)
    
    return  MVHM


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

    
    forceCoupledFD[nodes-1] = 0
    forceCoupledFD[nodes] = 0 
    forceCoupledFD[2*nodes-1] = 0
    forceCoupledFD[2*nodes] = 0

    uFDMVHM = solve(Coupling(nodes,h,x),forceCoupled)
    uFD = solve(CouplingFDFD(nodes,h,xFD),forceCoupledFD)
    uFDFull = solve(FDM(nodesFull,h),forceFull(nodesFull,h))
    
    uSlice = np.array(np.concatenate((uFDMVHM[0:nodes],uFDMVHM[nodes+3:2*nodes+2],uFDMVHM[2*nodes+5:len(x)])))
    uSliceFD = np.array(np.concatenate((uFD[0:nodes],uFD[nodes+1:2*nodes],uFD[2*nodes+1:len(x)])))

    uVHM = solve(VHM(nodesFull,h,xFull),forceFull(nodesFull,h))

    #print(max(uSlice),max(uSliceFD),max(uFDFull),max(uVHM))
    print(max(abs(uSlice-uSliceFD)))
    #print(np.array(np.concatenate((xFD[0:nodes],xFD[nodes+1:2*nodes],xFD[2*nodes+1:len(xFD)])))) 

    plt.axvline(x=1,c="#536872")
    plt.axvline(x=2,c="#536872")
    
    if example == "Quartic" or example == "Linear-cubic" or example =="Linear" or example == "Cubic" or example == "Quadratic":
        
        plt.plot(xFull,uSlice,label=r"$\delta$=1/"+str(int(n/2))+"",c="black",marker=markers[i-4],markevery=n)
        plt.plot(xFull,uSliceFD,label=r"$\delta$=1/"+str(int(n/2))+"",c="red",marker=markers[i-8],markevery=n)
        plt.ylabel("Error in displacement w.r.t. FDM")
        plt.ylabel("Error in displacement w.r.t. FDM")

        #plt.plot(xFull,uFDFull,label=r"$\delta$=1/"+str(int(n/2))+"",c="blue",marker=markers[i-4],markevery=n)
        #plt.plot(xFull,uVHM,label=r"$\delta$=1/"+str(int(n/2))+"",c="blue",marker=markers[i-4],markevery=n)
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
    file = "con_neumann_"+ str(c) + ".csv"
    # print(file)
    np.savetxt(file, con, delimiter=",")


