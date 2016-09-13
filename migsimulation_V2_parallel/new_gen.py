#! /usr/bin/env python
## This program creates a text data file
## which contains the location of point scatters.
#####################################################
# Created by Xin Liu on Nov 14, 2010
# Indiana University
# This program is under GNU public license
#####################################################

import os
import sys
import math

def rotate(x, y, rotang):
	# usage:
	# xp,yp = rotate(x, y, rotang)
	# rotation around z axis!
	x1=math.cos(rotang)*x-math.sin(rotang)*y
	y1=math.sin(rotang)*x+math.cos(rotang)*y
	return x1, y1
def sph2cart(theta, phi, r):
	#% theta: longitude
	#% phi: latitude
	#% r: radius (length)
	X=math.cos(phi)*math.cos(theta)
	Y=math.cos(phi)*math.sin(theta)
	Z=math.sin(phi)
	return X,Y,Z

def cart2sph(X, Y, Z):
	#% theta: longitude
	#% phi: latitude
	#% r: radius (length)
	if X>0:
	    theta=math.atan(Y/X)
	else:
	    theta=math.atan(Y/X)+math.pi
	
	r=math.sqrt(X*X+Y*Y+Z*Z)
	phi=math.atan(Z/math.sqrt(X*X+Y*Y))
	
	#[theta, phi, r]=cart2sph( X , Y  , Z)
	#theta=r2d(theta)
	#phi=r2d(phi)
	return theta,phi,r

def coordrot(phi, theta,lat0, lon0, xAz):
	# all angles in radians.
	#% phi: latitude
	#% theta: longitude
	# lat0, lon0: the coords of the origin points.
	# xAz: the azimuth of 
	x1,x2,x3=sph2cart(theta-lon0, phi, 1)
	alpha=math.pi/2-xAz
	beta=lat0
	#alpha=0
	#beta=0 ##debug only
	#print 'alpha=', alpha, 'beta=', beta
	x2p,x3p=rotate(x2, x3, alpha)
	x1p=x1
	x1pp,x3pp=rotate(x1p, x3p, beta)
	x2pp=x2p;
	lonreal,latreal,r=cart2sph(x1pp, x2pp, x3pp)	
	lonreal=math.fmod(lonreal+lon0, 2*math.pi)
	#print 'old', phi, theta
	#print 'new',latreal,lonreal 
	return latreal, lonreal

	
#outputfl="pointsourcelocs.dat"
outputfl="pointsourcelocs.tmp"
try:
	fsout=open(outputfl, "w")
except IOError:
	print "The file does not exist, exiting gracefully"

if len(sys.argv) == 6:
	print "happy!"
	lat1=math.radians(float(sys.argv[1]))
	lon1=math.radians(float(sys.argv[2]))
	lat2=math.radians(float(sys.argv[3]))
	lon2=math.radians(float(sys.argv[4]))
	zp=float(sys.argv[5]);
else:
	#print "Error! Please type the command like this: "
	print 'test, test!'
	print "python generate_point_scatters.py lat lon latp lonp depth\n"
	#sys.exit(0)
# getting started!

dx=35.0; # the step in km 
R=6370; # Earth's radius.
ddeg=dx/(R); # convert to degree.
print dx, ddeg 
x1offset=20
x2offset=43
n1=60
n2=70
xAz=75.0
xAz=math.radians(xAz)
lon0=-117.4497
lat0=41.4230
lon0=math.radians(lon0)
lat0=math.radians(lat0)
zp=410
lat1=0-(x2offset)*ddeg
lat2=0+(n2-(x2offset+1))*ddeg
lon1=lon0-(x1offset)*ddeg
lon2=lon0+(n1-(x1offset+1))*ddeg

amp=1
lat=lat1
while lat <= lat2:
	lon=lon1
	while lon <= lon2:
		#line=[str(lat), str(lon), str(zp), str(amp)]
		#line= "%(la)f %(lo)f %(dep)f %(amp)f\n" % {'la': math.degrees(lat), 'lo': math.degrees(lon), 'dep':zp, 'amp':amp} 
		latnew,lonnew=coordrot(lat, lon,lat0, lon0, xAz)
		print math.degrees(latnew), math.degrees(lonnew)
		line= "%(la)f %(lo)f %(dep)f %(amp)f\n" % {'la': math.degrees(latnew), 'lo': math.degrees(lonnew), 'dep':zp, 'amp':amp} 
		fsout.writelines(line)
		lon=lon+ddeg
	lat=lat+ddeg

fsout.close()
		
