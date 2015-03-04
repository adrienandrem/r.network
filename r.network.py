#!/usr/bin/env python
#
############################################################################
#
# MODULE:      r.network
# CREATED:     2014-11-07T10:40:53-0300
# AUTHOR:      Adrien Andre - adr dot andre at laposte dot net
# PURPOSE:     Generate a vector network from a raster map
#
#############################################################################

#%Module
#% description: Creates network from elevation raster map
#% keyword: raster
#% keyword: elevation
#% keyword: network
#% keyword: vector
#%end
#%option G_OPT_R_INPUT
#%end
#%option G_OPT_V_OUTPUT
#%end
#%option
#% key: threshold
#% type: double
#% description: Maximum slope threshold
#% required : no
#% answer : 1.
#%end

import os
import sys

import grass.script as grass
from grass.exceptions import CalledModuleError

from grass.pygrass.gis.region import Region
from grass.pygrass.raster import RasterSegment
from grass.pygrass.vector import VectorTopo
from grass.pygrass.vector.geometry import Line

from math import sqrt

PROCESS = 1

usermask = None
mapset = None
pid = None
overwrite = True
quiet = False

class Vertex(object):

    def __init__(self, i, j, previous):
        self.i = i
        self.j = j
        self.previous = previous

# Current cell to its neighbors connection
class Arc(object):

    def __init__(self, source, target, slope): # , side_slope):
        nodes = [source, target]
        nodes.sort()
        self.source = nodes[0]
        self.target = nodes[1]
        self.slope      = slope
        #self.side_slope = side_slope

    def __repr__(self):
        return "Edge(%s, %s)" % (self.source, self.target)

    def __eq__(self, other):
        if isinstance(other, Arc):
            return ((self.source == other.source) and (self.target == other.target))
        else:
            return False
    def __ne__(self, other):
        return (not self.__eq__(other))

    def __hash__(self):
        return hash(self.__repr__())

    def __lt__(self, other):
        if (self.source == other.source):
            return (self.target < other.target)
        else:
            return (self.source < other.source)

    def __gt__(self, other) :
        if (self.source == other.source):
            return (self.target > other.target)
        else:
            return (self.source > other.source)


def build_arc(i, j, width, index, i_shift, j_shift, slope_forw_r): # , slope_side_r):
        slope_forw = abs(slope_forw_r[i, j])
        #slope_side = abs(slope_side_r[i, j])
        index_target = (i + i_shift)*width + j + j_shift

        return Arc(index, index_target, slope_forw) # , slope_side)


def network(elev, network, threshold):
    '''See: Epstein, 2006, https://dx.doi.org/10.1287/opre.1060.0331'''

    region = Region()
    x_res = region.ewres
    y_res = region.nsres

    if grass.find_file('MASK', mapset = mapset)['file']:
        usermask = "usermask_mask." + pid
        print "A user raster mask (MASK) is present. Saving it..."
        grass.run_command('g.rename', quiet = quiet, raster = ('MASK', usermask))
    else:
        print "No user raster mask (MASK) found."


    print "Computing slope."

    print "E"
    grass.mapcalc("slope_forw_e  = ({0}[ 0,  1] - {0}[0, 0])/ewres()".format(elev), overwrite = True, quiet = quiet)
    print "N"
    grass.mapcalc("slope_forw_n  = ({0}[ 1,  0] - {0}[0, 0])/nsres()".format(elev), overwrite = True, quiet = quiet)
    print "W"
    grass.mapcalc("slope_forw_w  = ({0}[ 0, -1] - {0}[0, 0])/ewres()".format(elev), overwrite = True, quiet = quiet)
    print "S"
    grass.mapcalc("slope_forw_s  = ({0}[-1,  0] - {0}[0, 0])/nsres()".format(elev), overwrite = True, quiet = quiet)

    print "NE"
    grass.mapcalc("slope_forw_ne = ({0}[-1,  1] - {0}[0, 0])/sqrt(ewres()*ewres() + nsres()*nsres())".format(elev), overwrite = True, quiet = quiet)
    print "NW"
    grass.mapcalc("slope_forw_nw = ({0}[-1, -1] - {0}[0, 0])/sqrt(ewres()*ewres() + nsres()*nsres())".format(elev), overwrite = True, quiet = quiet)
    print "SW"
    grass.mapcalc("slope_forw_sw = ({0}[ 1, -1] - {0}[0, 0])/sqrt(ewres()*ewres() + nsres()*nsres())".format(elev), overwrite = True, quiet = quiet)
    print "SE"
    grass.mapcalc("slope_forw_se = ({0}[ 1,  1] - {0}[0, 0])/sqrt(ewres()*ewres() + nsres()*nsres())".format(elev), overwrite = True, quiet = quiet)

    d = sqrt(y_res*y_res + 2.0*x_res*2.0*x_res) # FIXME: Wrong if ewres() != nsres()
    print "ENE"
    grass.mapcalc("slope_forw_ene = ({0}[-1,  2] - {0}[0, 0])/{1}".format(elev, d), overwrite = True, quiet = quiet)
    print "NNE"
    grass.mapcalc("slope_forw_nne = ({0}[-2,  1] - {0}[0, 0])/{1}".format(elev, d), overwrite = True, quiet = quiet)
    print "NNW"
    grass.mapcalc("slope_forw_nnw = ({0}[-2, -1] - {0}[0, 0])/{1}".format(elev, d), overwrite = True, quiet = quiet)
    print "WNW"
    grass.mapcalc("slope_forw_wnw = ({0}[-1, -2] - {0}[0, 0])/{1}".format(elev, d), overwrite = True, quiet = quiet)
    print "WSW"
    grass.mapcalc("slope_forw_wsw = ({0}[ 1, -2] - {0}[0, 0])/{1}".format(elev, d), overwrite = True, quiet = quiet)
    print "SSW"
    grass.mapcalc("slope_forw_ssw = ({0}[ 2, -1] - {0}[0, 0])/{1}".format(elev, d), overwrite = True, quiet = quiet)
    print "SSE"
    grass.mapcalc("slope_forw_sse = ({0}[ 2,  1] - {0}[0, 0])/{1}".format(elev, d), overwrite = True, quiet = quiet)
    print "ESE"
    grass.mapcalc("slope_forw_ese = ({0}[ 1,  2] - {0}[0, 0])/{1}".format(elev, d), overwrite = True, quiet = quiet)

#    print "Computing side slope."
#    print "East"
#    grass.mapcalc("slope_side_e  = (({0}[ 1,  0] + {0}[ 1,  1])/2.0 - ({0}[-1,  0] + {0}[-1,  1])/2.0)/(2.0*nsres())".format(elev), overwrite = True, quiet = quiet)
#    print "North"
#    grass.mapcalc("slope_side_n  = (({0}[ 0,  1] + {0}[-1,  1])/2.0 - ({0}[ 0, -1] + {0}[-1, -1])/2.0)/(2.0*ewres())".format(elev), overwrite = True, quiet = quiet)
#    print "West"
#    grass.mapcalc("slope_side_w  = (({0}[ 0, -1] + {0}[-1, -1])/2.0 - ({0}[ 1,  0] + {0}[ 1, -1])/2.0)/(2.0*nsres())".format(elev), overwrite = True, quiet = quiet)
#    print "South"
#    grass.mapcalc("slope_side_s  = (({0}[ 0, -1] + {0}[ 1, -1])/2.0 - ({0}[ 0,  1] + {0}[ 1,  1])/2.0)/(2.0*ewres())".format(elev), overwrite = True, quiet = quiet)
#    print "N-E"
#    grass.mapcalc("slope_side_ne = ({0}[ 0,  1] - {0}[-1,  0])/sqrt(ewres()*ewres() + nsres()*nsres())".format(elev), overwrite = True, quiet = quiet)
#    print "N-W"
#    grass.mapcalc("slope_side_nw = ({0}[-1,  0] - {0}[ 0, -1])/sqrt(ewres()*ewres() + nsres()*nsres())".format(elev), overwrite = True, quiet = quiet)
#    print "S-W"
#    grass.mapcalc("slope_side_sw = ({0}[ 0, -1] - {0}[ 1,  0])/sqrt(ewres()*ewres() + nsres()*nsres())".format(elev), overwrite = True, quiet = quiet)
#    print "S-E"
#    grass.mapcalc("slope_side_se = ({0}[ 1,  0] - {0}[ 0,  1])/sqrt(ewres()*ewres() + nsres()*nsres())".format(elev), overwrite = True, quiet = quiet)

    # restoring user's mask, if present
    if usermask:
        print "Restoring user mask (MASK)..."
        try:
            grass.run_command('g.rename', quiet=quiet, raster = (usermask, 'MASK'))
        except CalledModuleError:
            print "Failed to restore user MASK!"
        usermask = None


    region = Region() # g.region align=dem

    print "Fetching base maps..."
    dem   = RasterSegment(elev)
    mask  = RasterSegment("MASK")

    print "Reading base maps..."
    dem.open('r')
    mask.open('r')


    print "Fetching slope maps..."
    slope_forw_e  = RasterSegment('slope_forw_e')
    slope_forw_n  = RasterSegment('slope_forw_n')
    slope_forw_w  = RasterSegment('slope_forw_w')
    slope_forw_s  = RasterSegment('slope_forw_s')
    slope_forw_ne = RasterSegment('slope_forw_ne')
    slope_forw_nw = RasterSegment('slope_forw_nw')
    slope_forw_sw = RasterSegment('slope_forw_sw')
    slope_forw_se = RasterSegment('slope_forw_se')
    slope_forw_ene = RasterSegment('slope_forw_ene')
    slope_forw_nne = RasterSegment('slope_forw_nne')
    slope_forw_nnw = RasterSegment('slope_forw_nnw')
    slope_forw_wnw = RasterSegment('slope_forw_wnw')
    slope_forw_wsw = RasterSegment('slope_forw_wsw')
    slope_forw_ssw = RasterSegment('slope_forw_ssw')
    slope_forw_sse = RasterSegment('slope_forw_sse')
    slope_forw_ese = RasterSegment('slope_forw_ese')


    print "Reading slope maps..."
    slope_forw_e.open('r')
    slope_forw_n.open('r')
    slope_forw_w.open('r')
    slope_forw_s.open('r')
    slope_forw_ne.open('r')
    slope_forw_nw.open('r')
    slope_forw_sw.open('r')
    slope_forw_se.open('r')
    slope_forw_ene.open('r')
    slope_forw_nne.open('r')
    slope_forw_nnw.open('r')
    slope_forw_wnw.open('r')
    slope_forw_wsw.open('r')
    slope_forw_ssw.open('r')
    slope_forw_sse.open('r')
    slope_forw_ese.open('r')

#    print "Fetching side slope maps..."
#    slope_side_e  = RasterSegment('slope_side_e')
#    slope_side_n  = RasterSegment('slope_side_n')
#    slope_side_w  = RasterSegment('slope_side_w')
#    slope_side_s  = RasterSegment('slope_side_s')
#    slope_side_ne = RasterSegment('slope_side_ne')
#    slope_side_nw = RasterSegment('slope_side_nw')
#    slope_side_sw = RasterSegment('slope_side_sw')
#    slope_side_se = RasterSegment('slope_side_se')
#
#
#    print "Reading side slope maps..."
#    slope_side_e.open('r')
#    slope_side_n.open('r')
#    slope_side_w.open('r')
#    slope_side_s.open('r')
#    slope_side_ne.open('r')
#    slope_side_nw.open('r')
#    slope_side_sw.open('r')
#    slope_side_se.open('r')


    print "Building arcs..."
    # TODO: Build a cost function (power or exponential)
    arc_set = set()
    rows = region.rows
    cols = region.cols
    width = cols
    for i in range(2, rows - 2):
        for j in range(2, cols - 2):

            if mask[i, j] != PROCESS:
                continue

            index_source = i*width + j

            # TODO: Improve these 16 ugly similar calls. DRY!

            # East
            i_shift =  0
            j_shift =  1
            if mask[i + i_shift, j + j_shift] == PROCESS:
                arc_set.add(build_arc(i, j, width, index_source, i_shift, j_shift, slope_forw_e)) # , slope_side_e))

            # North
            i_shift = -1
            j_shift =  0
            if mask[i + i_shift , j + j_shift] == PROCESS:
                arc_set.add(build_arc(i, j, width, index_source, i_shift, j_shift, slope_forw_n)) # , slope_side_n))

            # West
            i_shift =  0
            j_shift = -1
            if mask[i + i_shift, j + j_shift] == PROCESS:
                arc_set.add(build_arc(i, j, width, index_source, i_shift, j_shift, slope_forw_w)) # , slope_side_w))

            # South
            i_shift =  1
            j_shift =  0
            if mask[i + i_shift, j + j_shift] == PROCESS:
                arc_set.add(build_arc(i, j, width, index_source, i_shift, j_shift, slope_forw_s)) # , slope_side_s))


            # North-East
            i_shift = -1
            j_shift =  1
            if mask[i + i_shift, j + j_shift] == PROCESS:
                arc_set.add(build_arc(i, j, width, index_source, i_shift, j_shift, slope_forw_ne)) # , slope_side_ne))

            # North-West
            i_shift = -1
            j_shift = -1
            if mask[i + i_shift, j + j_shift] == PROCESS:
                arc_set.add(build_arc(i, j, width, index_source, i_shift, j_shift, slope_forw_nw)) # , slope_side_nw))

            # South-West
            i_shift =  1
            j_shift = -1
            if mask[i + i_shift, j + j_shift] == PROCESS:
                arc_set.add(build_arc(i, j, width, index_source, i_shift, j_shift, slope_forw_sw)) # , slope_side_sw))

            # South-East
            i_shift =  1
            j_shift =  1
            if mask[i + i_shift, j + j_shift] == PROCESS:
                arc_set.add(build_arc(i, j, width, index_source, i_shift, j_shift, slope_forw_se)) # , slope_side_se))


            # ENE
            i_shift = -1
            j_shift =  2
            if mask[i + i_shift, j + j_shift] == PROCESS:
                arc_set.add(build_arc(i, j, width, index_source, i_shift, j_shift, slope_forw_ene))

            # NNE
            i_shift = -2
            j_shift =  1
            if mask[i + i_shift, j + j_shift] == PROCESS:
                arc_set.add(build_arc(i, j, width, index_source, i_shift, j_shift, slope_forw_nne))

            # NNW
            i_shift = -2
            j_shift = -1
            if mask[i + i_shift, j + j_shift] == PROCESS:
                arc_set.add(build_arc(i, j, width, index_source, i_shift, j_shift, slope_forw_nnw))

            # WNW
            i_shift = -1
            j_shift = -2
            if mask[i + i_shift, j + j_shift] == PROCESS:
                arc_set.add(build_arc(i, j, width, index_source, i_shift, j_shift, slope_forw_wnw))

            # WSW
            i_shift =  1
            j_shift = -2
            if mask[i + i_shift, j + j_shift] == PROCESS:
                arc_set.add(build_arc(i, j, width, index_source, i_shift, j_shift, slope_forw_wsw))

            # SSW
            i_shift =  2
            j_shift = -1
            if mask[i + i_shift, j + j_shift] == PROCESS:
                arc_set.add(build_arc(i, j, width, index_source, i_shift, j_shift, slope_forw_ssw))

            # SSE
            i_shift =  2
            j_shift =  1
            if mask[i + i_shift, j + j_shift] == PROCESS:
                arc_set.add(build_arc(i, j, width, index_source, i_shift, j_shift, slope_forw_sse))

            # ESE
            i_shift =  1
            j_shift =  2
            if mask[i + i_shift, j + j_shift] == PROCESS:
                arc_set.add(build_arc(i, j, width, index_source, i_shift, j_shift, slope_forw_ese))


    print "Built {0} arcs.".format(len(arc_set))

    print "Filtering ({0} threshold)...".format(threshold)
    arcs = [arc for arc in list(arc_set) if arc.slope < threshold]
    print "{0} arcs remaining.".format(len(arcs))

    print "Sorting..."
    arcs.sort()


    # http://grass.osgeo.org/grass71/manuals/libpython/pygrass_vector.html#working-with-vector-objects
    net = VectorTopo(network)
    columns = [(u"cat", "INTEGER PRIMARY KEY"),
               (u"source", "INTEGER"),
               (u"target", "INTEGER"),
               (u"slope", "DOUBLE PRECISION")] # ,
               # (u"side_slope", "DOUBLE PRECISION")]
    net.open('w', tab_name = network, tab_cols = columns)

    x_min = region.west
    y_max = region.north
    width = region.cols

    print "Building geometries..."
    for a in arcs:
        s_j = a.source % width
        s_i = (a.source - s_j)/width

        t_j = a.target % width
        t_i = (a.target - t_j)/width

        s_x = x_min + (0.5 + s_j)*x_res
        s_y = y_max - (0.5 + s_i)*y_res
        s_z = dem[s_i, s_j]

        t_x = x_min + (0.5 + t_j)*x_res
        t_y = y_max - (0.5 + t_i)*y_res
        t_z = dem[t_i, t_j]

        net.write(Line([(s_x, s_y, s_z), (t_x, t_y, t_z)]), (a.source, a.target, a.slope, )) # e.side_slope, ))

    net.table.conn.commit()
    net.close()


def main():
    '''Main function.

    r.network in=dem out=net'''

    # Workspace UID
    global usermask, mapset, pid, overwrite, quiet
    mapset = grass.gisenv()['MAPSET']
    pid = str(os.getpid())

    input  = options['input']
    output = options['output']
    threshold = float(options['threshold'])

    network(input, output, threshold)


if __name__ == '__main__':
    options, flags = grass.parser()
    sys.exit(main())
