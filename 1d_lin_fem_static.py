# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 07:02:33 2019

@author: Loukit
"""
import numpy as np


#Defining the parameters of the problem

#Geometric properties of the rod
Length = 5.
Area = 1.

#Elastic Material properties
mu = 50. #Shear Modulus
nu = 0.3 #Poisson's Ratio
const = 2*mu*Area*(1-nu)/(1-2*nu)

#Loading Details
bodyForce = 10.
traction = 2.

#Mesh Details
nEls = 10 #Number of elements
elType = 3 #Number of nodes per element (2 for linear element, 3 for quadratic element)
NNodes = (elType -1) * nEls +1 #Total number of nodes in the model


def assemble_mesh_data(nEls,elType,NNodes,Length):
    nodes = np.zeros(NNodes) #Nodal Coordinates
    el_connectivity = np.zeros((nEls,elType),dtype = int) #element connectivity matrix
    
    for i in range(0,NNodes):
        nodes[i] = Length*(i)/(NNodes-1)
    
    for i in range(0,nEls):
        if (elType ==3):
            el_connectivity[i,0] = 2*i
            el_connectivity[i,2] = 2*(i+1)-1
            el_connectivity[i,1] = 2*(i+1)
        elif (elType == 2):
            el_connectivity[i,0] = i+1
            el_connectivity[i,1] = i+2
            
    return nodes,el_connectivity

nodes,el_connectivity = assemble_mesh_data(nEls,elType,NNodes,Length)

#Integration points and weights for 2-point Gaussian quadrature
nPoints = elType-1
if (nPoints ==2):
    w = np.array([1,1],dtype = float)
    xi = np.array([-0.5773502692,0.5773502692])
elif (nPoints == 1):
    w = np.array([2,0],dtype = float)
    xi = np.array([0.,0.])


#Assembling global stiffness and force vectors
K = np.zeros((NNodes,NNodes))
f = np.zeros((NNodes,1))

for i in range(0,nEls):
    #Getting the coorindates of each node on the current element
    nodesCoords = np.zeros(elType)
    for j in range(0,elType):
        nodesCoords[j] = nodes[el_connectivity[i,j]]
    
    #Defining the element stiffness matrix and load vector for current element
    kEl = np.zeros([elType,elType])
    fEl = np.zeros(elType)
    
    for j in range(0,nPoints):
        #Computing N and dN/dxi at the current integration point
        N = np.zeros(elType)
        dNdxi = np.zeros(elType)
        
        if (elType == 3):
            N[0] = -0.5 * xi[j] * (1.0 - xi[j])
            N[1] = 0.5 * xi[j] * (1.0 + xi[j])
            N[2] = (1.0 - (xi[j]**2))
            dNdxi[0] = -0.5 + xi[j]
            dNdxi[1] = 0.5  + xi[j]
            dNdxi[2] = -2.0 * xi[j]
            
        elif (elType == 2):
            N[0] = 0.5 * (1.0 - xi[j])
            N[1] = 0.5 * (1.0 + xi[j])
            dNdxi[0] = -0.5
            dNdxi[1] = 0.5
        #Computing dx/dxi, J and dN/dx at the current integration point
        dxdxi = 0.
        for k in range(0,elType):
            dxdxi = dxdxi + dNdxi[k]*nodesCoords[k]
        J = np.abs(dxdxi)
        dNdx = np.zeros(elType)
        for k in range(0,elType):
            dNdx[k] = dNdxi[k]/dxdxi
        #Adding the contribution to the element stiffness and force vector from current integration point
        for k in range(0,elType):
            fEl[k] = fEl[k] + w[j]*bodyForce*J*N[k]
            for l in range(0,elType):
                kEl[k,l] = kEl[k,l] + const*w[j]*J*dNdx[k]*dNdx[l]
    
    #Adding the stiffness and residual from the current element into global matrices
    for p in range(0,elType):
        row = el_connectivity[i,p]
        f[row] = f[row]+ fEl[p]
        for q in range(0,elType):
            col = el_connectivity[i,q]
            K[row,col] = K[row,col] + kEl[p,q]
#Adding the extra forcing term from the traction at x=L
f[NNodes-1] = f[NNodes-1] + traction

#Enforcing Boundary Conditions u=0 at first node
for i in range(0,NNodes):
    K[0,i] = 0.
K[0,0] = 1.
f[0] = 0.

#Solving the equation
u = np.linalg.solve(K,f)

#Plotting the results
import matplotlib.pyplot as plt

plt.plot(nodes,u)

    



