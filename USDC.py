#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 11 18:56:44 2020

@author: ehealy
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from PIL import Image, ImageFont, ImageDraw

class trace_segment:

    rho = 260 # in nm*Ohm
    def __init__(self, L, d, w):
        """
        Creates trace segment object.
        L : length of trace segment in mm
        d : depth of trace segment in nm
        w : width of trace segment in um
        A : xsection area of trace in nm*um
        R : total resistance of segment in kOhm
        """
        self.L = L
        self.d = d
        self.w = w
        self.A = d*w
        self.R = self.L*self.rho/self.A
    
  
class trace_database:

    def __init__(self, fname):
        """
        Opens .csv database file 
        """
        self.data = np.genfromtxt(fname, delimiter = ',', names = True,
            dtype = None, encoding = None)

class route:

    def __init__(self, tdb, route_id):
        """
        Creates route objects.
        tdb : trace database object
        route_id : electrical line id
        """
        
        rho = 260 
        
        self.trace_segments = []
        idxs = tdb.data['Route'] == route_id
        for row in tdb.data[idxs]:
            L = float(row['Length'])
            d = float(row['Depth'])
            w = float(row['Width'])
            self.trace_segments.append(trace_segment(L, d, w))

        self.integrate_Ls()
        self.integrate_Rs()
        self.R = self.integrated_R[-1]
        self.name = route_id

    def integrate_Ls(self):
        out = np.zeros(len(self.trace_segments) + 1)
        running_sum = 0
        for i, ts in enumerate(self.trace_segments):
            running_sum += ts.L
            out[i + 1] = running_sum
        self.integrated_L = out
        return out
    
    def integrate_Rs(self):
        out = np.zeros(len(self.trace_segments) + 1)
        running_sum = 0
        for i, ts in enumerate(self.trace_segments):
            running_sum += ts.R
            out[i + 1] = running_sum
        self.integrated_R = out
        return out         
    
    def plot_RvsL(self):
        R = self.integrated_R
        L = self.integrated_L
        f = interpolate.interp1d(L, R)
        xnew = np.arange(0, L[-1])
        ynew = f(xnew)
        plt.plot(L, R, 'o', xnew, ynew, color='blue', linewidth=2, markersize=6)
        plt.xlabel('Length (mm)')
        plt.ylabel('Resistance (kOhm)')
        plt.title('R vs. L')
        plt.show()
        return

    
def get_shorts_between_bias_lines(R_a, R_b, R_c, R_d, R_e):
    M = [[1,1,0,0,0], [0,0,1,1,0], [1,0,1,0,1], [1,0,0,1,1], [0,1,1,0,1]]
    M = np.array(M)
    I = np.linalg.inv(M) 
    R = [R_a, R_b, R_c, R_d, R_e]
    return np.matmul(I, R)


def get_shorts_to_ground(R_cont, R_p, R_n):
    M = np.array([[1,1,0], [1,0,1], [0,1,1]])
    I = np.linalg.inv(M) 
    R = [R_cont, R_p, R_n]
    return np.matmul(I, R)

    


def triangulate_BL_short(r1, r2, R_a, R_b, R_c, R_d, R_e):
    R1 = r1.integrated_R
    R2 = r2.integrated_R
    L1 = r1.integrated_L
    L2 = r2.integrated_L
    A = str(r1.name)
    B = str(r2.name)
    f, g = interpolate.interp1d(R1, L1), interpolate.interp1d(R2, L2)
    h = get_shorts_between_bias_lines(R_a, R_b, R_c, R_d, R_e)
    m, n = R1[-1]-h[1], R2[-1]-h[3]
    xnew1 = np.arange(0, R1[-1], 0.1)
    ynew1 = f(xnew1)
    plt.rcParams['figure.dpi'] = 120
    plt.ylabel('Length (mm)')
    plt.xlabel('Resistance (kOhm)')
    plt.title('R vs. L')
    plt.plot(xnew1, ynew1, color='blue', linewidth=1, label=A)
    plt.plot(R1, L1, 'o', color='blue', markersize=3)
    plt.plot(h[0], f(h[0]), '*b', markersize=12)
    plt.plot(m, f(m), '*b', markersize=12)
    xnew2 = np.arange(0, R2[-1], 0.1)
    ynew2 = g(xnew2)
    plt.ylabel('Length (mm)')
    plt.xlabel('Resistance (kOhm)')
    plt.title('R vs. L')
    plt.plot(xnew2, ynew2, color='red', linewidth=1, label=B) 
    plt.plot(R2, L2, 'o', color='red', markersize=3) 
    plt.plot(h[2], g(h[2]), '*r', markersize=12)
    plt.plot(n, g(n), '*r', markersize=12)
    xnew3 = [m, n]
    ynew3 = [f(m), g(n)]
    xnew4 = [h[0], h[2]]
    ynew4 = [f(h[0]), g(h[2])]
    plt.plot(xnew3, ynew3, linestyle='dashed', color='black')
    plt.plot(xnew4, ynew4, linestyle='dashed', color='black')
    plt.legend()
    plt.show()
    shortA = float(f(h[0]))
    shortB = float((g(h[2])))
    print('Short on {} is near L = '.format(A) + str(round(shortA)))
    print('Short on {} is near L = '.format(B) + str(round(shortB)))
    return 


def triangulate_ground_short(r1, R_a, R_b, R_c):
    R1 = r1.integrated_R
    L1 = r1.integrated_L
    f = interpolate.interp1d(R1, L1)
    h = get_shorts_to_ground(R_a, R_b, R_c)
    m = R1[-1]-h[1]
    xnew1 = np.arange(0, R1[-1], 0.1)
    ynew1 = f(xnew1)
    plt.ylabel('Length (mm)')
    plt.xlabel('Resistance (kOhm)')
    plt.title('R vs. L')
    plt.plot(R1, L1, 'o', xnew1, ynew1, color='blue', linewidth=2, markersize=6)
    plt.plot(h[0], f(h[0]), 'oc')
    plt.plot(m, f(m), 'oc')
    plt.show()
    print('Short is likely near L = ' + str(np.round(f(h[0]))) + ' and ' + str(np.round(f(h[1]))))
    return 



def circuit_diagram(R1, R2, R3, R4, R5):
    resistances = [R1, R2, R3, R4, R5]
    img = Image.open("DC_shorts.png")
    draw = ImageDraw.Draw(img)
    my_font = ImageFont.truetype("/System/Library/Fonts/Supplemental/Arial.ttf", 50)
    caption = ImageFont.truetype("/System/Library/Fonts/Supplemental/Arial.ttf", 20)
    draw.text((50,50),"resistances in kOhms", font=caption, fill=(255,0,0))
    draw.text((400,80),str(round(resistances[0],1)), font=my_font, fill=(255,0,0))
    draw.text((700,80),str(round(resistances[1],1)), font=my_font, fill=(255,0,0))
    draw.text((400,530),str(round(resistances[2],1)), font=my_font, fill=(255,0,0))
    draw.text((700,530),str(round(resistances[3],1)), font=my_font, fill=(255,0,0))
    draw.text((730,300),str(round(resistances[4],1)), font=my_font, fill=(255,0,0))
    img.save("shorts.png")
    return 





# def triangulate_short2(r1, r2, R_a, R_b, R_c, R_d, R_e):
#     R1 = r1.integrated_R
#     R2 = r2.integrated_R
#     L1 = r1.integrated_L
#     L2 = r2.integrated_L
#     f, g = interpolate.interp1d(R1, L1), interpolate.interp1d(R2, L2)
#     h = get_shorts_between_bias_lines(R_a, R_b, R_c, R_d, R_e)
#     m, n = R1[-1]-h[1], R2[-1]-h[3]
#     xnew1 = np.arange(0, R1[-1], 0.1)
#     ynew1 = f(xnew1)
#     plt.ylabel('Length (mm)')
#     plt.xlabel('Resistance (kOhm)')
#     plt.title('R vs. L')
#     plt.plot(R1, L1, 'o', xnew1, ynew1, color='blue', linewidth=2, markersize=6)
#     xnew2 = np.arange(0, R2[-1], 0.1)
#     ynew2 = g(xnew2)
#     plt.plot(R2, L2, 'o', xnew2, ynew2, color='red', linewidth=2, markersize=6) 
#     xnew3 = [h[0],h[2]]
#     ynew3 = [f(h[0]), g(h[2])]
#     plt.plot(xnew3, ynew3, color='black', linestyle='dashed')
#     plt.plot(m, f(m), 'ok')
#     plt.plot(h[2], g(h[2]), 'ok')
#     plt.plot(n, g(n), 'ok')
#     plt.show()
#     return





