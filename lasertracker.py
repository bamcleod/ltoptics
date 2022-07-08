import numpy as np
import matplotlib.pyplot as pl
from sklearn import preprocessing
import pandas

#-------------------------------------------------------------------------------
# FIT CIRCLE 2D
# - Find center [xc, yc] and radius r of circle fitting to set of 2D points
# - Optionally specify weights for points
#
# - Implicit circle function:
#   (x-xc)^2 + (y-yc)^2 = r^2
#   (2*xc)*x + (2*yc)*y + (r^2-xc^2-yc^2) = x^2+y^2
#   c[0]*x + c[1]*y + c[2] = x^2+y^2
#
# - Solution by method of least squares:
#   A*c = b, c' = argmin(||A*c - b||^2)
#   A = [x y 1], b = [x^2+y^2]
#-------------------------------------------------------------------------------
def fit_circle_2d(x, y, w=[]):
    
    A = np.array([x, y, np.ones(len(x))]).T
    b = x**2 + y**2

    # Modify A,b for weighted least squares
    if len(w) == len(x):
        W = np.diag(w)
        A = np.dot(W,A)
        b = np.dot(W,b)
    
    # Solve by method of least squares
    c ,resid, rank, s = np.linalg.lstsq(A,b,rcond=None)

#    print(c,resid,rank,s)
    
    # Get circle parameters from solution c
    xc = c[0]/2
    yc = c[1]/2
    r = np.sqrt(c[2] + xc**2 + yc**2)
    return xc, yc, r


#-------------------------------------------------------------------------------
# RODRIGUES ROTATION
# - Rotate given points based on a starting and ending vector
# - Axis k and angle of rotation theta given by vectors n0,n1
#   P_rot = P*cos(theta) + (k x P)*sin(theta) + k*<k,P>*(1-cos(theta))
#-------------------------------------------------------------------------------
def rodrigues_rot(P, n0, n1):
    
    # If P is only 1d array (coords of single point), fix it to be matrix
    if P.ndim == 1:
        P = P[np.newaxis,:]
    
    # Get vector of rotation k and angle theta
    n0 = n0/np.linalg.norm(n0)
    n1 = n1/np.linalg.norm(n1)
    k = np.cross(n0,n1)
    k = k/np.linalg.norm(k)
    theta = np.arccos(np.dot(n0,n1))
    
    # Compute rotated points
    P_rot = np.zeros((len(P),3))
    for i in range(len(P)):
        P_rot[i] = P[i]*np.cos(theta) + np.cross(k,P[i])*np.sin(theta) + k*np.dot(k,P[i])*(1-np.cos(theta))

    return P_rot


######
# Given 3 points compute a set of perpendicular vectors of length 1:
# i points from v1 to v2
# j is normal to the plane defined by v1, v2, v3
# k is perpendicular to i and j
#
def normalvectors(x, y):
    i = x
    k = np.cross(x,y)
    i = preprocessing.normalize(i, norm='l2')
    k = preprocessing.normalize(k, norm='l2')
    j = np.cross(k,i)
    return i, j , k

    
#######
# Transform a set of points to a new coordinate system
# pointlist is [x, y, z] or [[x1, y1, z1], ..., [xn, yn, zn]
# origin is [x0, y0, z0] the origin of the new coordinate system, expressed as measred in old CS
# i, j, k are the unit vectors of the new coordinate system, expressed as measured in the old CS
#
def transform(pointlist, origin, i, j, k, verbose=False):

    if verbose:
        print('Transforming', pointlist)
        print('with origin', origin)
        print('and vectors', i,j,k)
        
    new = []
    newx = []
    newy = []
    newz = []
    
    onepoint = False
    
    if len(pointlist.shape) == 1:
        onepoint = True
        pointlist = [pointlist]
    
    pointarray = np.array(pointlist)

    for row in range(len(pointlist)):
        if len(pointarray.shape) == len(origin.shape):
            tmp = pointlist[row] - origin[row]
        else:
            tmp = pointlist[row] - origin

        if (len(i.shape) == len(pointarray.shape)):
            ii = i[row]
            jj = j[row]
            kk = k[row]
        else:
            ii = i
            jj = j
            kk = k

        new.append([np.dot(tmp,ii), np.dot(tmp,jj), np.dot(tmp,kk)])
        newx.append(np.dot(tmp, ii) )
        newy.append(np.dot(tmp, jj) )
        newz.append(np.dot(tmp, kk) )

    new = np.array(new)
    newx = np.array(newx)
    newy = np.array(newy)
    newz = np.array(newz)
    
    result = np.transpose(np.array([newx, newy, newz]))
#    print(new)
    if onepoint:
#        print('onepoint', new[0])
        return (new[0])

    
    return new

"""
baselist is a list of indices into the list of retros defining which measurements define the baseplate coordinate system
"""
class LaserTrackerData:
    def __init__(self, filename, number_of_retros, extras=[], swaps=[], baselist=[1,0,2], xaxisindices=[], yaxisindices=[], origin=0):

        # Header contains degree symbol represented as 0xb0 so we need this encoding
        csvdata = pandas.read_csv(filename, sep=';', decimal=',', encoding='ISO-8859-1')

        X, Y, Z = 'X  [mm]', 'Y  [mm]', 'Z  [mm]'
        
        # Put into numpy array
        points = np.transpose(np.array([csvdata[X], csvdata[Y], csvdata[Z]]))

        # Replace any bad measurements with replacements taken at the end
        for extra in extras:
            points[extra[0]] = points[extra[1]]

        for swap in swaps:
            points[swap[0]], points[swap[1]] = points[swap[1]], points[swap[0]].copy()

        # Truncate to complete sets of data
        number_of_sets = int(len(points) / number_of_retros)
        points = points[:number_of_sets * number_of_retros]

        self.raw = np.transpose(points.reshape((number_of_sets, number_of_retros,3)), [1,0,2])

        # points raw[i][j][k]: retro i; data set j; k=x,y,z coords  in laser tracker coordinates

        #
        # Base plate wobble must be subtracted
        #

        # For each stage position transform the stage position into the coordinate system defined by the base plate.
        # The midpoint of positions 0 and 1 is the origin
        # 1 - 0 defines the X axis
        # normal to the baseplate plane is the Z axis
        # perpendicular to those two is the Y axis
        
        i, j, k = normalvectors(self.raw[xaxisindices[1]]-self.raw[xaxisindices[0]], 
                                self.raw[yaxisindices[1]]-self.raw[yaxisindices[0]])

        self.normal = [i,j,k]

        # Now transform the probe coordinates into the base plate coords

        self.new = []
        for point in self.raw:
            self.new.append(transform(point, self.raw[origin], i, j, k))

        self.new = np.array(self.new)        

######################
# Find the center of the circle defined by the stage retro
# Do this in 3D because the bearing plane is not exactly parallel to the baseplate
# Use the algorithm described here
# https://meshlogic.github.io/posts/jupyter/curve-fitting/fitting-a-circle-to-cluster-of-3d-points/
#####################

# Returns:
#  xc, yc: the center of the circle
#  r: the radius of the circle
#  normal:  normal vector
#  ang: the angles of the points projected onto the circle plane
#  rad: the radii of the points projected onto the circle plane
#
def CircleFit (points):
    
    #-------------------------------------------------------------------------------
    # (1) Fitting plane by SVD for the mean-centered data
    # Eq. of plane is <p,n> + d = 0, where p is a point on plane and n is normal vector
    #-------------------------------------------------------------------------------
    P_mean = points.mean(axis=0)
    P_centered = points - P_mean
    U,s,V = np.linalg.svd(P_centered)
    # Normal vector of fitting plane is given by 3rd column in V
    # Note linalg.svd returns V^T, so we need to select 3rd row from V^T
    normal = V[2,:]
    d = -np.dot(P_mean, normal)  # d = -<p,n>

    #-------------------------------------------------------------------------------
    # (2) Project points to coords X-Y in 2D plane
    #-------------------------------------------------------------------------------
    P_xy = rodrigues_rot(P_centered, normal, [0,0,1])

    # Find the center of the circle defined by the stage retro
    #

    xc, yc, r = fit_circle_2d(P_xy[:,0], P_xy[:,1])
    print ("Centered points: ", P_centered)
    print ("Projected points: ", P_xy)
    print ("Center of circle:", xc, yc)
    print ("Radius of circle:", r)

    # Compute the angle of each point, measured from -y axis
    #
    ang = np.degrees(np.arctan2(P_xy[:,0] - xc, -(P_xy[:,1] - yc)))
    rad = np.sqrt((P_xy[:,0]-xc)**2 + (P_xy[:,1]-yc)**2)

    # Translate back to original coordinate system

    circle_center = rodrigues_rot(np.array([xc, yc, 0]), [0,0,1], normal) + P_mean

    return circle_center[0], r, normal, ang, rad, P_mean


#######################
# Compute rigid body transformation between 2 sets of 3D points
# 
# Description: https://nghiaho.com/?page_id=671
# Code: http://nghiaho.com/uploads/code/rigid_transform_3D.py_

# Input: expects Nx3 matrix of points
# Returns R,t
# R = 3x3 rotation matrix
# t = 3x1 column vector

def rigid_transform_3D(A, B):

    print (A)
    print (B)
    assert len(A) == len(B)

    N = A.shape[0]; # total points

    centroid_A = np.mean(A, axis=0)
    centroid_B = np.mean(B, axis=0)
    
    # centre the points
    AA = A - np.tile(centroid_A, (N, 1))
    BB = B - np.tile(centroid_B, (N, 1))

    # dot is matrix multiplication for array
    H = np.transpose(AA) @ BB

    U, S, Vt = np.linalg.svd(H)

    R = Vt.T @ U.T

    # special reflection case
    if np.linalg.det(R) < 0:
        #       print ("Reflection detected")
        Vt[2,:] *= -1
        R = Vt.T @ U.T

    t = -R @ centroid_A.T + centroid_B.T

    #    print t

    return R, t.reshape((-1,1))


def makealltheplots(lt, msrconfig, refconfig, resultscalefactor, ang_cmd, labeloffsets, titlestr, rootname, figaxes=None):
    results = {}
    configs = lt.keys()

    # Select only the moving retros for further analysis
    for config in configs:
        results[config] = lt[config].new[3:]

    xbase,ybase = lt[refconfig].new[:3,0,:2].T

    # Compute 1/2 the difference between +/- gravity
    diff = (results[msrconfig] - results[refconfig]) / resultscalefactor


    if figaxes == None:
        fig,axes = pl.subplots(3,2,figsize=(8,10),sharex=True,sharey=True)
    else:
        fig,axes = figaxes
        
    pl.axis('equal')
    # Put them all on one plot
    coordlabel = ['X', 'Y', 'Z']
    ang_color = ['r', 'g', 'b']
    for coord in [0, 1, 2]:
        ax = axes[coord,0]
        ax.plot(xbase,ybase,'go')
        for ang_i in range(len(ang_cmd)):
            good = (abs(diff[:,ang_i,0]) < 5.) & (abs(diff[:,ang_i,1]) < 5.) & (abs(diff[:,ang_i,2]) < 5. )
            ax.plot(results[refconfig][:,ang_i,0][good], results[refconfig][:,ang_i,1][good],color = ang_color[ang_i], marker = '+', linestyle='-', label=str(ang_cmd[ang_i]))
            for j in range(len(diff)):
                x = results[refconfig][j,ang_i,0]
                y = results[refconfig][j,ang_i,1]
                delta = diff[j,ang_i,coord]
                if good[j]:
                    deltastr = (1000 * delta).astype(int).astype(str)
                    ax.annotate(deltastr,(x+30,y+labeloffsets[ang_i]),color=ang_color[ang_i])
        #ax.set_axis('equal')
        ax.set_xlim(-600,1000)
        ax.set_ylim(-665,655)
        #        ax.legend(loc='upper left')
        #ax.title(titlestr + coordlabel[coord] + " deflection")
        #        pl.savefig(rootname + '-deflect-' + coordlabel[coord]+'.png')


    for config in [msrconfig]:
        if config == refconfig:
            continue
        print("Refconfig = ", refconfig, "Config = ", config)
        for ang_i in range(len(ang_cmd)):
            good = (abs(diff[:,ang_i,0]) < 5.) & (abs(diff[:,ang_i,1]) < 5.) & (abs(diff[:,ang_i,2]) < 5. )
            #
            # Transform the five points on the hub measured at 180deg to match the same points measured at 0deg
            #
            hub = results[refconfig][:,0,0]<500
            
            R,t = rigid_transform_3D(results[refconfig][hub,ang_i,:][good[hub]], results[config][hub,ang_i,:][good[hub]])
            A2 = (R @ results[refconfig][:,ang_i,:].T) + np.tile(t, (1, len(results[refconfig])))
            A2 = A2.T
            for coord in [0, 1, 2]:
                ax = axes[coord,1]
                ax.plot(results[refconfig][:,ang_i,0][good], results[refconfig][:,ang_i,1][good], color = ang_color[ang_i], marker = '+', linestyle='-', label=str(ang_cmd[ang_i]))
                for j in range(len(A2)):
                    x = results[refconfig][j,ang_i,0]
                    y = results[refconfig][j,ang_i,1]
                    delta = (results[config][j,ang_i,coord] - A2[j,coord]) / resultscalefactor
                    if good[j]:
                        deltastr = (1000 * delta).astype(int).astype(str)
                        ax.annotate(deltastr,(x+30,y+labeloffsets[ang_i]),color=ang_color[ang_i])

            R = R / resultscalefactor
            t = t / resultscalefactor
            print ("\nAngle = " + str(ang_cmd[ang_i]),end = ' ')
            print ("ThetaX: ", (R[1,2] * 1e6).astype(int),end=' ')
            print ("ThetaY: ", (R[0,2] * 1e6).astype(int),end=' ')
            print ("ThetaZ: ", (R[0,1] * 1e6).astype(int),end=' ')

            print ("Translation [micron]", end=' ')
            print ("Tx ", (t[0] * 1e3).astype(int), end=' ')
            print ("Ty ", (t[1] * 1e3).astype(int), end=' ')
            print ("Tz ", (t[2] * 1e3).astype(int), end=' ')

            print ("Probetip [micron]", end=' ')
            print ('X ', ((results[config][~hub,ang_i,0] - A2[~hub,0]).mean() / resultscalefactor * 1000).astype(int), end=' ')
            print ('Y ', ((results[config][~hub,ang_i,1] - A2[~hub,1]).mean() / resultscalefactor * 1000).astype(int), end=' ')
            print ('Z ', ((results[config][~hub,ang_i,2] - A2[~hub,2]).mean() / resultscalefactor * 1000).astype(int), end=' ')
            print (' ')

            # Print them again with no labels, so we can paste into a table
            print (str(ang_cmd[ang_i]),end = ' ')
            print ((R[1,2] * 1e6).astype(int),end=' ')
            print ((R[0,2] * 1e6).astype(int),end=' ')
            print ((R[0,1] * 1e6).astype(int),end=' ')
            print ((t[0][0] * 1e3).astype(int), end=' ')
            print ((t[1][0] * 1e3).astype(int), end=' ')
            print ((t[2][0] * 1e3).astype(int), end=' ')
            # probe tip
            for coord in [0,1,2]:
                print (((results[config][~hub,ang_i,coord] - A2[~hub,coord]).mean() / resultscalefactor * 1000).astype(int), end=' ')
            print (" ")
            
    axes[0][0].legend(ncol=3,loc='upper center',bbox_to_anchor=(0.5,0))
    axes[0][0].set_title("Raw deflections")
    axes[0][1].set_title("Fit residuals")
    for coord in [0,1,2]:
        ax = axes[coord,1]
        ax.text(1.05, 0.5, coordlabel[coord] ,
                rotation=0, size=16,
                bbox=dict(edgecolor='black', facecolor='none', pad=10, linewidth=2),
                ha='left', va='center', transform=ax.transAxes)
    
    fig.suptitle(titlestr,size=16)
    fig.savefig(rootname + '.png')
    
