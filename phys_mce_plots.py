import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import coordinateTransform as Lxmap
#from mpl_toolkits.axes_grid1 import make_axes_locatable

#Define a triangular grid that is useful for 1 pixel
def translateTriangles(deltaX,deltaY):
    x = np.asarray([0,1,1,0,0.5]) + deltaX
    y = np.asarray([0,0,1,1,0.5]) + deltaY
    triangles = [[0,4,1],[3,4,0], [2,4,3], [1,4,2]]
    return x,y,triangles

#Define a triangular grid by replicating that one pixel
def triangleGrid(numX, numY):
    totalX = np.array([])
    totalY = np.array([])
    totalTriangle = np.empty((0,3), int)
    for indexX in np.arange(0,numX):
        for indexY in np.arange(0, -numY, -1):
            indexAdj = len(totalX)
            
            x,y,triangles = translateTriangles(indexX, indexY)
            totalX = np.append(totalX,x)
            totalY = np.append(totalY,y)
            totalTriangle = np.append( totalTriangle, np.add(triangles,indexAdj), 0)
            

    return totalX, totalY, totalTriangle


x,y,triangles = triangleGrid(4,4)



triang = mtri.Triangulation(x,y, triangles)
#fig,axs = plt.subplots(nrows = 1, ncols=1)
#axs = axs.flatten()

#zfaces is the modification if you want to plot something (like opt efficiency , etc.). 

#zfaces = np.arange(0,len(triangles))
zfaces = np.zeros(len(triangles))

def indexzface(R,C,band,pol):
    col = int(C)
    row = int(R)
    pixelloc = 16*(col - 1) + 4*(row - 1)
    if band == 'U' and pol == 'A':
        pixel = pixelloc + 2
    elif band == 'U' and pol == 'B':
        pixel = pixelloc + 3
    elif band == 'L' and pol == 'A':
        pixel = pixelloc
    else:
        pixel = pixelloc + 1
    return pixel

def plotPhysBowtie(mce_matrix, titlename = "", axis = []):
    for indexRow in np.arange(0,33):
        for indexCol in np.arange(0,2):
            #Get the phys detector location
            physPixel = bowtiemap.getMCE2Phys(indexRow, indexCol)
            physRow = physPixel[0]
            physCol = physPixel[1]
            physBand = physPixel[2]
            physPol = physPixel[3]
            #Convert to flatten index
            index = indexzface(physRow, physCol, physBand, physPol)
            #Assign to correct value
            zfaces[index] = mce_matrix[indexRow, indexCol]


    fig, ax = plt.subplots()
    plt.title(titlename)
    #cs = ax.tripcolor(x,y,triangles, facecolors = zfaces,edgecolors = 'k')

    if len(axis) > 0:
        cs = ax.tripcolor(x,y,triangles, facecolors = zfaces,edgecolors = 'k', vmin = axis[0], vmax = axis[1])
            
    else:

        cs = ax.tripcolor(x,y,triangles, facecolors = zfaces,edgecolors = 'k')
    cbar = fig.colorbar(cs, shrink=0.9, pad = 0.0)
    ax.axis('off')
    #Hide the frame
    #ax.spines["top"].set_visible(False)
    #ax.spines["right"].set_visible(False)
    #ax.spines["bottom"].set_visible(False)
    #ax.spines["left"].set_visible(False)
    ax.set_aspect(1)

    ax.text(2.0, -3.5, 'Det Col', horizontalalignment = 'center', verticalalignment = 'bottom', color ='k', fontsize = 12, fontweight = 'bold')
    ax.text(0.5, -3.2, '1', horizontalalignment = 'center', verticalalignment = 'bottom', color = 'k', fontsize = 10, fontweight = 'bold')
    ax.text(1.5, -3.2, '2', horizontalalignment = 'center', verticalalignment = 'bottom', color = 'k', fontsize = 10, fontweight = 'bold')
    ax.text(2.5, -3.2, '3', horizontalalignment = 'center', verticalalignment = 'bottom', color = 'k', fontsize = 10, fontweight = 'bold')
    ax.text(3.5, -3.2, '4', horizontalalignment = 'center', verticalalignment = 'bottom', color = 'k', fontsize = 10, fontweight = 'bold')


    ax.text(-0.5, -1.0, 'Det Row', horizontalalignment = 'center', verticalalignment = 'center', color = 'k', fontsize = 12, fontweight = 'bold', rotation = 90)
    ax.text(-0.2, 0.5, '1', horizontalalignment = 'left', verticalalignment = 'center', color = 'k', fontsize = 10, fontweight = 'bold')
    ax.text(-0.2, -0.5, '2', horizontalalignment = 'left', verticalalignment = 'center', color = 'k', fontsize = 10, fontweight = 'bold')
    ax.text(-0.2, -1.5, '3', horizontalalignment = 'left', verticalalignment = 'center', color = 'k', fontsize = 10, fontweight = 'bold')
    ax.text(-0.2, -2.5, '4', horizontalalignment = 'left', verticalalignment = 'center', color = 'k', fontsize = 10, fontweight = 'bold')

    
    ax.text(0.5, 0.33, 'LA', horizontalalignment = 'center', verticalalignment = 'center', color = 'k', fontsize = 7)
    
    ax.text(0.5, 0.66, 'UA', horizontalalignment = 'center', verticalalignment = 'center', color = 'k', fontsize = 7)

    
    ax.text(0.33, 0.5, 'LB', horizontalalignment = 'center', verticalalignment = 'center', color = 'k', fontsize = 7)

    
    ax.text(0.66, 0.5, 'UB', horizontalalignment = 'center', verticalalignment = 'center', color = 'k', fontsize = 7)

    #plt.savefig(filename)
    #plt.show()
