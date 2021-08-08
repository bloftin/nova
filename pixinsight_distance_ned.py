# The goal of this script is an aid to find the farthest objects in your astophotography pictures
# This script is to take a list of objects detected via annotation in pixinisight and dumped to a file 
# called All_objects.txt and it connects to NED online database to get redshift info on objects
# and then calcualte the Hubble distance based on redshift and then creates a custom catalogue for pixinsight so you can annotate your images with the distance
# of objects (Pixinisght_custom_catalogue.txt)
# Copyright (c) 2021 Benjamin Loftin.
# 
# This program is free software: you can redistribute it and/or modify  
# it under the terms of the GNU General Public License as published by  
# the Free Software Foundation, version 3.
#
# This program is distributed in the hope that it will be useful, but 
# WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License 
# along with this program. If not, see <http://www.gnu.org/licenses/>.

from astroquery.ned import Ned
from math import *
import astropy.units as u
import astropy.coordinates as coord
import astropy.units as u
from astropy.coordinates import SkyCoord, ICRS, Galactic

fout    = open("All_objects_distances.txt", "w")
fout_pi = open("Pixinsight_custom_catalogue.txt","w")

# write header
fout.write('Name, Hubble Distance (Mly), RA, DEC, Velocity CMB (km/s), Redshift (Helio)\n')
fout_pi.write('RA' + "\t" + 'DEC' + "\t" + "NAME" + "\t" + "DIAMETER"  + "\n")

with open('All_objects.txt', 'r') as reader:
    line = reader.readline()
    while line != '':  # The EOF char is an empty string
        print(line, end='')
        line = reader.readline()

        astro_object = ''

        if line:
            astro_object = line.split(';')[0]
            ra_pi_deg = line.split(';')[1]
            dec_pi_deg = line.split(';')[2]

        
        if  astro_object:

            # Check to see if PGC type, if so need to search with LEDA instad
            if "PGC" in astro_object:
                astro_object = astro_object.replace("PGC", "LEDA ")

            # ['photometry'|'positions'|'diameters'|'redshifts'|'references'|'object_notes'].
            # result_table = Ned.get_table("NGC5229", table='redshifts', verbose=True)
            try:
                result_table_pos   = Ned.get_table(astro_object, table='positions', verbose=True)
                result_table_red   = Ned.get_table(astro_object, table='redshifts', verbose=True)
                result_table_diam  = Ned.get_table(astro_object, table='diameters', verbose=True)


                # RA ranges from 0-24 hours due to earth rotation 24 hours
                # DEC ranges -90 to +90, +90 deg is the north celestial pol, 0 cel equator, -90 south
                # so need to deal with DEC having +/- in string
                ra_hms_str = result_table_pos['RA'][0].decode('ASCII')
                dec_hms_str = result_table_pos['DEC'][0].decode('ASCII')

                hour = float(ra_hms_str.split('h')[0])
                minute = float(ra_hms_str.split('h')[1].split('m')[0])
                second = float(ra_hms_str.split('h')[1].split('m')[1].split('s')[0])

                ra_deg = abs(hour) + abs(minute)/60 + abs(second/3600)
                #take care of sign
                ra_deg =  ra_deg*copysign(1, ra_deg)
                # multiply by earth rate, 1 hour = 15 deg
                ra_deg = ra_deg* 15.0

                # if greater than 180, subtract 360
                #if ra_deg > 180.0:
                #    ra_deg = ra_deg-360.0

                ra_rad = radians(ra_deg)

                dec = float(dec_hms_str.split('d')[0])
                minute = float(dec_hms_str.split('d')[1].split('m')[0])
                second = float(dec_hms_str.split('d')[1].split('m')[1].split('s')[0])

                dec_deg = abs(dec) + abs(minute)/60 + abs(second/3600)
                #take care of sign
                dec_deg =  dec_deg*copysign(1, dec_deg)
                dec_rad = radians(dec_deg)

                # Heliocentric redshift 
                z_helio = float(result_table_red['Published Redshift'][0])

                # Calcualte Diameter needed for pixinsight annotation, convet arcsec to arcmin
                # which is what we get back from NED, assuming diameter to encompass ellipse is 2*major_axis
                diam = 2.0*float(result_table_diam['NED Major Axis'][0])/60.0

                # convert to cosmic background radiation reference to get hubble distance
                # calculation is from NED docs, https://ned.ipac.caltech.edu/Documents/Guides/Calculators
                # Heliocentric to 3K Background 	264.14 deg 	+48.26 deg 	371.0 km/sec 	ApJ 473, 576, 1996
                I_apex = radians(264.14)
                b_apex = radians(48.26)
                V_apex = 371.0
                hubble_constant =  67.8    # km/s
                Mpc2ly = 3261563.7769443   # 1 Mpc = 3261563.7769443 ly 
                c = 299792.458   # speed of light in km/s
                vz_helio = c*z_helio     # multiply by speed of light to get velocity in helio frame

                # need to convert the helio (equatorial coordinates) to galactic for below velocity conversion
                # https://en.wikipedia.org/wiki/Galactic_coordinate_system
                icrs = SkyCoord(ra=ra_rad*u.radian, dec=dec_rad*u.radian, frame='icrs')
                gc = icrs.galactic

                #v_cmb = v_helio + V_apex*(sin(ra_rad)*sin(b_apex) + cos(ra_rad)*cos(b_apex)*cos(dec_rad-I_apex))
                v_cmb = vz_helio + V_apex*(sin(gc.b.radian)*sin(b_apex) + cos(gc.b.radian)*cos(b_apex)*cos(gc.l.radian-I_apex))

                hubble_distance_Mpc = v_cmb / hubble_constant
                hubble_distance_ly = hubble_distance_Mpc*Mpc2ly

                print("RA Helio (deg):       "  + str(ra_deg))
                print("DEC Helio (deg):       "  + str(dec_deg))
                print("l Galactic (deg):       "  + str(gc.l.degree))
                print("b Galactic (deg):       "  + str(gc.b.degree))
                print("V_helio (km/s):  " + str(vz_helio))
                print("V_CMB (km/s):    "  + str(v_cmb))
                print("Hubble Distance (Mpc): "  + str(hubble_distance_Mpc))
                print("Hubble Distance (ly): "  + str(hubble_distance_ly))
                print("Hubble Distance (Mly): "  + str(hubble_distance_ly/1E6))
                print("Diameter (arcsec): " + str(diam))

                #output data to file
                fout.write(astro_object + ', ' + str(hubble_distance_ly/1E6) + ', ' + ra_hms_str + ', ' + 
                           dec_hms_str + ', ' + str(v_cmb) + ', ' + str(z_helio) +'\n')

                print_dist = "{:.0f}".format(hubble_distance_ly/1E6)
                fout_pi.write( ra_pi_deg + "\t" + dec_pi_deg + "\t" + astro_object + " , " + print_dist + " Mly" + "\t" + "{:.2f}".format(diam)  + "\n")

            except:
                print("Data not avail: " + astro_object)