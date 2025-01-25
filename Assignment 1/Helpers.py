import numpy as np
from matplotlib import pyplot as plt
import math

### Helper function relevant to cylindrical coordinates.
def cylindricalToCartesian(R, phi, Z):
    """
    Conversion of cylindrical coordinates to Cartesian coordinates.
    """
    x = R * np.cos(phi)
    y = R * np.sin(phi)
    z = Z
    
    return x, y, z
    
def cartesianToCylindrical(x, y, z):
    """
    Conversion of Cartesian coordinates to cylindrical coordinates.
    """
    R = np.sqrt(x**2 + y**2)
    phi = np.arctan(y / x)
    Z = zfill
    
    return R, phi, Z

def cylindricalCross(R1, phi1, Z1, R2, phi2, Z2):
    """
    Cross product in cylindrical coordinates. https://www.quora.com/What-is-a-cross-product-in-cylindrical-coordinates-What-are-some-examples
    Is identical to normal cross product.
    """
    R = phi1 * Z2 - Z1 * phi2
    phi = -(R1 * Z2 - Z1 * R2)
    Z = R1 * phi2 - phi1 * R2
    
    return R, phi, Z

def cylindricalInner(R1, phi1, Z1, R2, phi2, Z2):
    """
    Inner product in cylindrical coordinates.
    """
    x1, y1, z1 = cylindricalToCartesian(R1, phi1, Z1)
    x2, y2, z2 = cylindricalToCartesian(R2, phi2, Z2)
    
    return x1 * x2 + y1 * y2 + z1 * z2

def cylindricalAddition(R1, phi1, Z1, R2, phi2, Z2):
    """
    Addition of cylindrical coordinates. https://math.stackexchange.com/questions/1365622/adding-two-polar-vectors
    """
    R = np.sqrt(R1**2 + R2**2 + 2 * R1 * R2 * np.cos(phi2 - phi1))
    phi = phi1 + np.arctan2(R2 * np.sin(phi2 - phi1), R1 + R2 * np.cos(phi2 - phi1))
    Z = Z1 + Z2
    
    return R, phi, Z

### Helper methods.
def BStar(BR, Bphi, BZ, R, q, dt, m, VphiCurrent, VphiOld):
    """
    Method for converting a normal magnetic field into a new magnetic field which uses the generalized toroidal momentum.
    """
    BStarR = q * dt * BR / m
    BStarphi = q * dt * Bphi / m
    BStarZ = q * dt * (BZ + m * (1.5 * VphiCurrent + 0.5 * VphiOld)/(q * R)) / m
    
    return BStarR, BStarphi, BStarZ

def gyrationTime(m, B, Z=2, e=1.602*(10**-19)):
    """
    https://en.wikipedia.org/wiki/Plasma_parameters and 3MF100 and Freidberg
    """
    return (m*2*np.pi)/(Z*e * B)

def makePlot(dataSets, R0=1, names=0, AXReturn=False):
    """
    Makes a R-Z and a X-Y plot of the particle trajectory.
    """
    if names==0:
        names = ["Particle " + str(i + 1) for i in range(len(dataSets))]
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize = (18,9))
    zAbsoluteMax = 0
    absoluteMax = 0
    for i in range(len(dataSets)):
        data = dataSets[i]
        amountOfDatapoints = len(data["R"])
        
        xData = []
        yData = []
        zData = []
        for j in range(amountOfDatapoints):
            x, y, z = cylindricalToCartesian(data["R"][j], data["phi"][j], data["Z"][j])
            xData.append(x)
            yData.append(y)
            zData.append(z)
        
        ax1.plot(data["R"], data["Z"], label=names[i])
        ax1.plot([data["R"][amountOfDatapoints - 1]], [data["Z"][amountOfDatapoints - 1]], "o")#, label="Particle")
        ax1.set_title("R-Z plane")
        ax1.set_xlabel("R (meter)")
        ax1.set_ylabel("Z (meter)")
        zAbsoluteMaxNew = max([abs(min(data["Z"])), abs(max(data["Z"]))])
        if zAbsoluteMax < zAbsoluteMaxNew:
            zAbsoluteMax = zAbsoluteMaxNew
        ax1.set_ylim([-1.1*zAbsoluteMax, 1.1*zAbsoluteMax])
        
        #fig, ax = plt.subplots(1, 1, figsize = (8,8))
        ax2.plot(xData, yData, label=names[i])
        ax2.plot(xData[amountOfDatapoints - 1], yData[amountOfDatapoints - 1], "o")#, label="Particle")
        ax2.set_title("X-Y plane")
        ax2.set_xlabel("X (meter)")
        ax2.set_ylabel("Y (meter)")
        
        xAbsoluteMax = max([abs(min(xData)), abs(max(xData))])
        yAbsoluteMax = max([abs(min(yData)), abs(max(yData))])
        absoluteMaxNew = 1.1*max([xAbsoluteMax, yAbsoluteMax])
        if absoluteMax< absoluteMaxNew:
            absoluteMax = absoluteMaxNew
        ax2.set_xlim([-absoluteMax, absoluteMax])
        ax2.set_ylim([-absoluteMax, absoluteMax])
    
    ax1.plot([R0], [0], "o", label="Magnetic axis")
    ax2.plot(np.sin(np.linspace(0, 2*np.pi, num=100)), np.cos(np.linspace(0, 2*np.pi, num=100)), label="Magnetic axis")
    ax2.plot([0], [0], "o")
    ax1.legend()
    ax2.legend()
    if AXReturn:
        return fig, ax1, ax2
    else:
        plt.show()
    

### Fields that could be used.
def BField3(R, phi, Z, R0, B0=5, Bp0=1):
    """
    The magnetic field in a tokamak as defined in assignment 1 question 3.
    """
    r = np.sqrt((R - R0)**2 + Z**2)
    
    BR = Bp0 * Z * R0 / (r * R)
    Bphi = B0 * R0 / R
    BZ = Bp0 * (R0 - R) * R0 / (r * R)
    
    return BR, Bphi, BZ

def EField3(R, phi, Z, R0, E0=0, Ep0=0):
    """
    The electric field in a tokamak as defined in assignment 1 question 3.
    """
    #r = np.sqrt((R - R0)**2 + Z**2)
    
    return 0, 0, 0

def EField5(R, phi, Z, R0, E0=5, Ep0=1):
    """
    The electric field in a tokamak as defined in assignment 1 question 5.
    """
    r = np.sqrt((R - R0)**2 + Z**2)
    
    ER = E0 * (R - R0) / r
    Ephi = 0
    EZ = E0 * Z / r
    
    return ER, Ephi, EZ