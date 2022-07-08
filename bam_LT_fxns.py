from lasertracker import *

import numpy as np
import time



def Transform(origin_of_new_coord_sys, x_direction_vector_of_new_coord_sys, y_direction_vector_of_new_coord_sys, z_direction_vector_of_new_coord_sys):

    def trans(points):
        return transform(points, origin_of_new_coord_sys, x_direction_vector_of_new_coord_sys, y_direction_vector_of_new_coord_sys, z_direction_vector_of_new_coord_sys)

    origin_of_old_coord_sys = trans(np.array([0,0,0]))
    x_direction_vector_of_old_coord_sys = trans(np.array([1,0,0])) - origin_of_old_coord_sys
    y_direction_vector_of_old_coord_sys = trans(np.array([0,1,0])) - origin_of_old_coord_sys
    z_direction_vector_of_old_coord_sys = trans(np.array([0,0,1])) - origin_of_old_coord_sys

    def invtrans(points):
        return transform(points, origin_of_old_coord_sys, x_direction_vector_of_old_coord_sys, y_direction_vector_of_old_coord_sys, z_direction_vector_of_old_coord_sys)
    
    return trans, invtrans


class measurement_set:
    def __init__(self, LTcoords_bench_origin_smr, LTcoords_bench_xaxis_smr, LTcoords_bench_yaxis_smr, Zemax_origin_in_benchcoordinates, Zemax_unitvectors_in_bench_coordinates, zemax_prescrip_fname='OneFlat_oapdesign_using35-596_prescription_REFSURFZERO.txt'):

        fileObject = open(zemax_prescrip_fname, "r")
        data = fileObject.read()
        self.data = data
        self.zemax_coordinates = ZemaxData(data)
        normals_zemax = self.zemax_coordinates.normal
      
        bench_origin_inLTcoords, lb_i, lb_j, lb_k, LT_origin_inbenchcoords, LT_i, LT_j, LT_k = establish_bench_coordinates(LTcoords_bench_xaxis_smr,LTcoords_bench_origin_smr,LTcoords_bench_yaxis_smr)

        self.transform_LT_to_bench,    self.transform_bench_to_LT    = Transform(bench_origin_inLTcoords, lb_i, lb_j, lb_k)

        self.transform_bench_to_Zemax, self.transform_Zemax_to_bench = Transform(Zemax_origin_in_benchcoordinates,
                                                                       Zemax_unitvectors_in_bench_coordinates[0],
                                                                       Zemax_unitvectors_in_bench_coordinates[1],
                                                                       Zemax_unitvectors_in_bench_coordinates[2])

        self.normals_bench = self.transform_Zemax_to_bench(normals_zemax) - self.transform_Zemax_to_bench(np.array([0,0,0]))

        self.bench_origin_in_Zemaxcoordinates = self.transform_bench_to_Zemax(np.array([0,0,0]))

        self.bench_pointlist = self.transform_Zemax_to_bench(self.zemax_coordinates.vertex)
        self.bench_origin_inLTcoords = bench_origin_inLTcoords
        self.lb_i = lb_i
        self.lb_j = lb_j
        self.lb_k = lb_k
        self.LT_origin_inbenchcoords = LT_origin_inbenchcoords
        self.LT_i = LT_i
        self.LT_j = LT_j
        self.LT_k = LT_k
        self.Zemax_origin_in_benchcoordinates = Zemax_origin_in_benchcoordinates



    def select_surface(self, surface_index, use_back_srf=False):
        surfnum = int(surface_index)
        self.surface_coords_bench = self.bench_pointlist[surfnum]
        if not use_back_srf:
            self.surface_normal = self.normals_bench[surfnum]
        else:
            self.surface_normal = -1*self.normals_bench[surfnum]


        self.transform_Zemax_to_optic, self.transform_optic_to_Zemax = Transform(self.zemax_coordinates.vertex[surfnum],
                                                                                 self.zemax_coordinates.rotation[surfnum].T[0],
                                                                                 self.zemax_coordinates.rotation[surfnum].T[1],
                                                                                 self.zemax_coordinates.rotation[surfnum].T[2])
            
        print('Selected surface', surface_index, 'With normal', self.surface_normal)
    
    def measure_SMR_real(self, msrd_SMR_pos_LTcoord, cartesian=True):
        if cartesian:
            msrd_SMR_pos_LTcoord_XYZ = msrd_SMR_pos_LTcoord
        else: # convert from HVD to XYZ first
            msrd_SMR_pos_LTcoord_XYZ = LT_HVD_to_XYZ(msrd_SMR_pos_LTcoord[0], msrd_SMR_pos_LTcoord[1], msrd_SMR_pos_LTcoord[2])

        msrd_SMR_pos_benchcoord = transform(msrd_SMR_pos_LTcoord_XYZ, self.bench_origin_inLTcoords, self.lb_i, self.lb_j, self.lb_k)
        
        self.measured_SMR_position_benchcoord = msrd_SMR_pos_benchcoord
    
    def get_virtualimg_posn(self, cartesian=True):
        dist_bench = self.measured_SMR_position_benchcoord - self.surface_coords_bench 
        d_N = np.dot(dist_bench,self.surface_normal)
        print('dist dot normal', d_N)
        self.virtual_SMR_position_benchcoord = self.measured_SMR_position_benchcoord - 2*d_N*self.surface_normal
        self.virtual_SMR_posn_LTcoord_Cart = transform(self.virtual_SMR_position_benchcoord, self.LT_origin_inbenchcoords, self.LT_i, self.LT_j, self.LT_k)
        self.virtual_SMR_posn_LTcoord_sph = LT_XYZ_to_HVD(self.virtual_SMR_posn_LTcoord_Cart[0], self.virtual_SMR_posn_LTcoord_Cart[1], self.virtual_SMR_posn_LTcoord_Cart[2])
        if cartesian:
            print('Virtual image LT coords, XYZ:', self.virtual_SMR_posn_LTcoord_Cart)
            
        else:
            print('Virtual image LT coords, HVD:', self.virtual_SMR_posn_LTcoord_sph)
     
    def compare_VImsrmt_VIfiducial(self, VI_measured_posn, cartesian=True):
        if cartesian:
            d = VI_measured_posn - self.virtual_SMR_posn_LTcoord_Cart
           
        else:
            VI_msr_Cart= LT_HVD_to_XYZ(VI_measured_posn[0],VI_measured_posn[1],VI_measured_posn[2])
            d = VI_msr_Cart - self.virtual_SMR_posn_LTcoord_Cart
            
        self.last_meas_distance = np.linalg.norm(d)
        self.last_meas_distvector = d
        print('Distance (x,y,z [LT coords])', d)
        print('3D distance [mm]', np.linalg.norm(d))                                           
                                                       
    def get_measured_mirrornormal(self, VI_measured_posn, cartesian=True):
        if cartesian:
            VI_LT_xyz = VI_measured_posn
        else:
            VI_LT_xyz = LT_HVD_to_XYZ(VI_measured_posn[0], VI_measured_posn[1], VI_measured_posn[2])
            
        vi_msred_benchcoord = transform(VI_LT_xyz, self.bench_origin_inLTcoords, self.lb_i, self.lb_j, self.lb_k)
        VI_to_real = vi_msred_benchcoord - self.measured_SMR_position_benchcoord
        measurednorm = VI_to_real/np.linalg.norm(VI_to_real)
        print('Measured normal, bench coordinates', measurednorm)
        dotpr = np.dot(measurednorm, self.surface_normal)
        theta = np.sqrt(1 - (dotpr**2))
        print('Angle between:', theta)

def LT_HVD_to_XYZ(h,v,d):
    v_rad = v*np.pi/180
    h_rad = h*np.pi/180
    y = d*np.sin(v_rad)*np.cos(h_rad)
    x = d*np.sin(v_rad)*np.sin(h_rad)
    z = d*np.cos(v_rad)
    return np.array([x,y,z])

def LT_XYZ_to_HVD(x,y,z):
    d = np.sqrt(x**2 + y**2 + z**2)
    v = np.arccos(z/d)
    if x < 0 and y < 0:
        h = np.arctan(y/x)
        h_deg = -1*(180-(-1*(h*180/np.pi)-90))
    else:
        h = np.arctan(-y/x)
        h_deg = -1*(180-(-1*(-1*(h*180/np.pi)-90)))
        
    if np.abs(h_deg) > 180:
        if h_deg < 0:
            h_deg = h_deg+180
        else:
            h_deg = h_deg-180
    v_deg = v*180/np.pi
    return np.array([h_deg, v_deg, d])

def establish_bench_coordinates(smr1_LTcoord_xyz, smr2_LTcoord_xyz, smr3_LTcoord_xyz):
    '''
    use 3 measured SMR positions from LT to calculate origins + unit vectors for conversions from bench <-> LT coords
    INPUTS: [3]
    smr1_LTcoord_xyz, smr2_LTcoord_xyz, smr3_LTcoord_xyz
        Measured positions of 3 smrs in laser tracker coordinates
        
    OUTPUT: [8 outputs]
    bench_origin_inLTcoords, lb_i, lb_j, lb_k, LT_origin_inbenchcoords, LT_i, LT_j, LT_k
        bench_origin_inLTcoords: location of the bench origin in LT coordinates
        LT_origin_inbenchcoords: location of the laser tracker origin in Bench coordinates
        lb_i, lb_j, lb_k: the unit vectors used for LT -> Bench transformation
        LT_i, LT_j, LT_k: the unit vectors used for Bench -> LT transformation
    '''
    smr_LTcoords_xyz= np.array([smr1_LTcoord_xyz, smr2_LTcoord_xyz, smr3_LTcoord_xyz])

    bench_origin_inLTcoords, lb_i, lb_j, lb_k = PrepforTransform_LaserTrackerXYZ_to_BenchXYZ(smr_LTcoords_xyz)
    smr_Benchcoords = transform(smr_LTcoords_xyz, bench_origin_inLTcoords, lb_i, lb_j, lb_k)
#     print('SMR Positions, in bench coordinates:', new_smr_coords)

    # Prep for conversions between bench <-> Laser Tracker coords
    LT_x = np.array([1,0,0])
    LT_y = np.array([0,1,0])
    LT_z = np.cross(LT_y, LT_x) / np.linalg.norm(np.cross(LT_y, LT_x))

    # Conversion to bench coords part 1
    LT_i_pt = transform(LT_x, bench_origin_inLTcoords, lb_i, lb_j, lb_k)
    LT_j_pt = transform(LT_y, bench_origin_inLTcoords, lb_i, lb_j, lb_k)
    LT_k_pt = transform(LT_z, bench_origin_inLTcoords, lb_i, lb_j, lb_k)

    # Location of LT origin in Bench coordinates
    LT_origin_inbenchcoords = transform(np.array([0,0,0]), bench_origin_inLTcoords, lb_i, lb_j, lb_k)

    # Laser tracker unit vectors - correctly defined so they transform the SMR posns back to their original LT coords
    LT_i = (LT_i_pt - LT_origin_inbenchcoords)
    LT_j = (LT_j_pt-LT_origin_inbenchcoords)
    LT_k = (LT_origin_inbenchcoords - LT_k_pt)
    
    return bench_origin_inLTcoords, lb_i, lb_j, lb_k, LT_origin_inbenchcoords, LT_i, LT_j, LT_k

class ZemaxData:
    def __init__(self, readData):
        
        '''
        Input: readData, a Zemax .txt output file (prescription data) that has been read in with python's builtin open and read fxns
        nSurf: number of surfaces in the model
        beginning_of_data: beginning row; the one that says "Surf     R11     R12      R13     X     Tilt X"

        (The normal vector is -1*[R13, R23, R33])
        '''

        all_lines = readData.split('\n')

        # Find the global coordinate data
        #
        import re

        p = re.compile('^Surfaces')
        for line in all_lines:
            if p.match(line):
                number_of_surfaces = int(line.split(':')[1])
                print ('Nuber of surfaces', number_of_surfaces)
                nSurf = number_of_surfaces + 1
                break

        p = re.compile('^GLOBAL.*COORDINATES')
        for i, line in enumerate(all_lines):
            if p.match(line):
                beginning_of_data = i + 8
                print('Beginning of global coords:', beginning_of_data)

        lines = all_lines[beginning_of_data:beginning_of_data+(5*(nSurf-1))] # 5*nsurf-1

        self.comment  = [None] * nSurf
        self.vertex   = np.zeros((nSurf,3))
        self.rotation = np.zeros((nSurf,3,3))
        self.normal   = np.zeros((nSurf,3))
        self.nSurf    = nSurf

        for i, line in enumerate(lines):
            s1 = lines[i].split('\t')
            if len(s1) == 7:
                s2 = lines[i+1].split('\t')
                s3 = lines[i+2].split('\t')
                snum = int(s1[0])
                self.comment[snum] = s1[5]
                R11 = float(s1[1])
                R12 = float(s1[2])
                R13 = float(s1[3])
                X   = float(s1[4])
                R21 = float(s2[1])
                R22 = float(s2[2])
                R23 = float(s2[3])
                Y   = float(s2[4])
                R31 = float(s3[1])
                R32 = float(s3[2])
                R33 = float(s3[3])
                Z   = float(s3[4])
                self.vertex[snum] = np.array([X, Y, Z])
                self.rotation[snum] = np.array([[R11, R12, R13], [R21, R22, R23], [R31, R32, R33]])
                self.normal[snum] = np.array([-R13, -R23, -R33])

def PrepforTransform_LaserTrackerXYZ_to_BenchXYZ(smr_benchpositions_xyz, which_smr_origin=1):
    '''
    OUTPUT: origin, i, j, k
    '''
    smr1_benchpos_xyz = smr_benchpositions_xyz[0]
    smr2_benchpos_xyz = smr_benchpositions_xyz[1]
    smr3_benchpos_xyz = smr_benchpositions_xyz[2]
    
    origin = smr_benchpositions_xyz[int(which_smr_origin)] 
    diff12 = (smr1_benchpos_xyz - smr2_benchpos_xyz)
    diff32 = (smr3_benchpos_xyz - smr2_benchpos_xyz)
    i = diff12 / np.linalg.norm(diff12)
    k = np.cross(i, diff32) / np.linalg.norm(np.cross(i,diff32)) # Perpendicular to bench plane
    j = np.cross(k,i)
#     j = diff32 / np.linalg.norm(diff32)
#     k = np.cross(j,i) / np.linalg.norm(np.cross(j,i))
    return origin, i, j, k

def PrepforTransform_ZemaxXYZ_to_BenchXYZ(Zemax_origin_in_benchcoordinates, bench_i, bench_j, bench_k):
    '''
    output: bench_origin_in_Zemaxcoordinates
    '''
    bench_origin_in_Zemaxcoordinates = np.array(transform(np.array([0,0,0]), Zemax_origin_in_benchcoordinates, bench_k, bench_j, bench_i))
    return bench_origin_in_Zemaxcoordinates
