# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 15:04:00 2019

@author: 91917
"""
import os
from math import sqrt
import numpy as np
import matplotlib.pyplot as plt

DATA_PATH = os.path.join("input_files","fem_conststrain")


def read_data(data_path,file_name):
    try:
        file_path = os.path.join(data_path,file_name)
        file = open(file_path)
    except FileNotFoundError:
        print("File Not found. Please check the file name.")
        return None
    except:
        print("Unable to open file.")
        return None
        
    #Reading material properties
    try:
        mat_props = file.readline().strip().split()
        print(mat_props)
        E = float(mat_props[0])
        nu = float(mat_props[1])
        
        #Reading number of nodes and nodal coordinates
        NNodes = np.genfromtxt(file,dtype = int,max_rows = 1)
        NodeCoords = np.genfromtxt(file,max_rows = NNodes)
        
        #Reading number of elements and element connectivity matrix
        NEls = np.genfromtxt(file,dtype=int,max_rows=1)
        el_connect = np.genfromtxt(file,dtype = int,max_rows = NEls)
    
        #Reading number of prescribed displacements and prescribed displacements
        nfix = np.genfromtxt(file,dtype=int,max_rows=1)
        fixnodes = np.genfromtxt(file,max_rows = nfix)
        
        #Reading number of loaded elements and given loads
        ndload = np.genfromtxt(file,dtype=int,max_rows=1)
        if(ndload == 0):
            dloads = np.zeros((1,4))
        else:
            dloads = np.genfromtxt(file,max_rows = ndload).reshape((ndload,4))
        return E,nu,NNodes,NodeCoords,NEls,el_connect,nfix,fixnodes,ndload,dloads
#    except:
        print("File Reading error. Please check the format of input file.")
        return None
    

def element_stiffness_asm(xa,ya,xb,yb,xc,yc,Dmat):
    """
    Function to compute the element stiffness matrix for a constant strain
    triangle given the coordinates of the nodes a,b,c and D matrix 
    """

    area = (1/2)*abs( (xb-xa)*(yc-ya) - (xc-xa)*(yb-ya) )
    
    Bmat = compute_b_matrix(xa,ya,xb,yb,xc,yc)
    
    K_el = area*np.matmul(np.transpose(Bmat),np.matmul(Dmat,Bmat))
    return K_el

def el_force_vector(xa,ya,xb,yb,tx,ty):
    length = sqrt((xa-xb)**2 + (yb-ya)**2)
    rel = np.array([tx,ty,tx,ty])*0.5*length
    return rel
def compute_b_matrix(xa,ya,xb,yb,xc,yc):
    nax = -(yc-yb)/( (ya-yb)*(xc-xb) - (xa-xb)*(yc-yb) )
    nay =  (xc-xb)/( (ya-yb)*(xc-xb) - (xa-xb)*(yc-yb) )
    nbx = -(ya-yc)/( (yb-yc)*(xa-xc) - (xb-xc)*(ya-yc) )
    nby =  (xa-xc)/( (yb-yc)*(xa-xc) - (xb-xc)*(ya-yc) )
    ncx = -(yb-ya)/( (yc-ya)*(xb-xa) - (xc-xa)*(yb-ya) )
    ncy =  (xb-xa)/( (yc-ya)*(xb-xa) - (xc-xa)*(yb-ya) )
    
    Bmat = np.array([[nax,0,nbx,0,ncx,0],[0,nay,0,nby,0,ncy],
                     [nay,nax,nby,nbx,ncy,ncx]])
    
    return Bmat

def elastic_strain_mat(el_node_coords,u_el):
    xa = el_node_coords[0,0]
    ya = el_node_coords[0,1]
    xb = el_node_coords[1,0]
    yb = el_node_coords[1,1]
    xc = el_node_coords[2,0]
    yc = el_node_coords[2,1]
    uax = u_el[0,0]
    uay = u_el[0,1]
    ubx = u_el[1,0]
    uby = u_el[1,1]
    ucx = u_el[2,0]
    ucy = u_el[2,1]
    Bmat = compute_b_matrix(xa,ya,xb,yb,xc,yc)
    uel = np.array([uax,uay,ubx,uby,ucx,ucy])
    strain = np.matmul(Bmat,uel)
    return strain

def write_output(filename,u,NNodes,el_connect,Dmat):
    out_file = open(filename,'w')
    out_file.write("Nodal Displacements: \n")
    out_file.write("Node    u1       u2 \n")
    for i in range(0,NNodes):
        out_file.write('%d %8.4f %8.4f\n'% (i,u[2*i],u[2*i+1]))
    out_file.write('\n Strains and Stresses \n')
    out_file.write(' Element   e_11      e_22     e_12      s_11       s_22      s_12\n')
    smax = 0
    smises = np.zeros(NEls)
    for i in range(0,NEls):
        el_nodes = el_connect[i,:]
        el_node_coords = np.zeros((3,2))
        u_el = np.zeros((3,2))
        for j in range(0,3):
            el_node_coords[j,0] = NodeCoords[el_nodes[j],0]
            el_node_coords[j,1] = NodeCoords[el_nodes[j],1]
            u_el[j,0] = u[(2*el_nodes[j])]
            u_el[j,1] = u[(2*el_nodes[j])+1]
        strain = elastic_strain_mat(el_node_coords,u_el)
        stress = np.matmul(Dmat,strain)
        out_file.write('%d %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f\n'%(i,strain[0],strain[1],strain[2],stress[0],stress[1],stress[2]))
        smises[i] = np.sqrt((stress[0]-stress[1])**2 + 2*stress[2]**2)
        if(smises[i]>smax):
            smax = smises[i]
    out_file.close()


        
E,nu,NNodes,NodeCoords,NEls,el_connect,nfix,fixnodes,ndload,dloads = read_data(DATA_PATH,"FEM_conststrain_holeplate.txt")


#Plotting the undeformed mesh
# plt.figure(1)
# plt.triplot(NodeCoords[:,0],NodeCoords[:,1],el_connect)
# plt.show()

#Defining D matrix for computing element stiffness matrix
Dmat = np.array([[1-nu,nu,0],[nu,1-nu,0],[0,0,(1-2*nu)]])*(E/((1+nu)*(1-2*nu)))

#Assembling the global stiffness matrix
K_g = np.zeros((2*NNodes,2*NNodes),dtype = float) #Initializing

for el in range(0,NEls):
    a = el_connect[el,0]
    b = el_connect[el,1]
    c = el_connect[el,2]
    
    K_el = element_stiffness_asm(NodeCoords[a,0],NodeCoords[a,1],
                                 NodeCoords[b,0],NodeCoords[b,1],
                                 NodeCoords[c,0],NodeCoords[c,1],Dmat)
    
    for i in range(0,3):#Number of nodes per element
        for ii in range(0,2): #2d element. 2 dimension disp. vector
            for j in range(0,3):#number of nodes per element
               for jj in range(0,2):
                   row = 2*el_connect[el,i]+ii
                   col = 2*el_connect[el,j]+jj
                   K_g[row,col] = K_g[row,col] + K_el[(2*i+ii),(2*j+jj)]


#Assembling the global residual vector
resid = np.zeros(2*NNodes)
pointer = np.array([1,2,0])
for i in range(0,ndload):
    el = int(dloads[i,0])
    face = int(dloads[i,1])
    a = el_connect[el,face]
    b = el_connect[el,pointer[face]]
    r = el_force_vector(NodeCoords[a,0],NodeCoords[a,1],
                        NodeCoords[b,0],NodeCoords[b,1],
                        dloads[i,2],dloads[i,3])
    
    resid[2*a] += r[0]
    resid[(2*a)+1] += r[1]
    resid[2*b] += r[2]
    resid[(2*b)+1] += r[1]

#Including constraints in the stiffness matrix and residual force vector
for i in range(0,nfix):
    rw =int( (2*fixnodes[i,0])+fixnodes[i,1])
    for j in range(0,2*NNodes):
        K_g[rw,j]=0
    K_g[rw,rw] = 1.0
    """Ensuring that the corresponding u value becomes the value
    of the constraint. The equation effectively becomes
    K_g[rw,rw]*u[rw] = resid[rw], where K_g[rw,rw]=1.0(Rest of 
    the rows become zeros). Thus, u[rw] = resid[rw]"""
    resid[rw] = fixnodes[i,2]

#Solving for the displacement vector
u = np.linalg.solve(K_g,resid)

#computing the stress and strains and writing the output in a file
write_output("fem_constrain_holeplate_op.txt",u,NNodes,el_connect,Dmat)

#Plotting the results
plt.figure(1)
plt.triplot(NodeCoords[:,0],NodeCoords[:,1],el_connect)

defrmd_coords = NodeCoords
for i in range(0,NNodes):
    x = NodeCoords[i,0]
    y = NodeCoords[i,1]
    defrmd_coords[i,0] = x + u[2*i]
    defrmd_coords[i,1] = y + u[(2*i)+1]
    
plt.triplot(defrmd_coords[:,0],defrmd_coords[:,1],el_connect)
plt.show()


