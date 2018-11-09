#!/usr/bin/env python3
# I really hate python already:
from __future__ import print_function
import glob, os
import h5py
import argparse

class bcolors:
        HEADER = '\033[95m'
        OKBLUE = '\033[94m'
        OKGREEN = '\033[92m'
        WARNING = '\033[93m'
        FAIL = '\033[91m'
        ENDC = '\033[0m'
        BOLD = '\033[1m'
        UNDERLINE = '\033[4m'

def strictly_increasing(L):
    return all(x<y for x, y in zip(L, L[1:]))


def get_dset_name( fname ):
    from os.path import basename
    dset_name = basename(fname)
    dset_name = dset_name[0:dset_name.find('_')]

    return dset_name


def get_timestamp( fname ):
    import re
    from os.path import basename

    fname = basename(fname)
    # extract everything between "_" and "." so mask_00000.h5 gives 00000
    # note escape character "\"
    m = re.search('\_(.+?)\.', fname )

    if m:
        timestamp = m.group(1)
    else:
        print("An error occured: we couldn't extract the timestamp")

    return timestamp


def print_list( l ):
    for p in l:
        print(bcolors.HEADER + p + " " + bcolors.ENDC, end='')
    # print just one newline
    print(' ')


def warn( msg ):
    print( bcolors.FAIL + "WARNING! " + bcolors.ENDC + msg)


def uniquelist( l ):
    # if the list has only one unique value, return it, if not, error
    if len(l) == 0:
        return None
    l = sorted(list(set(l)))
    if len(l) == 1:
        return(l[0])
    elif len(l) == 0:
        warn('uniquelist: List is completely empty!')
        return None
    else:
        warn('uniquelist: List ist not unique...something went wrong.')
        print('these are the values we found in the list:')
        print(l)
        return l[0]


def write_xmf_file_wabbit(args, outfile, times, timestamps, prefixes, scalars, vectors, directory):
    print('-------------------------')
    print('- WABBIT module         -')
    print('-------------------------')

    fid = open(outfile, 'w')
    # use any file to get blocksize and dimensionality
    file = directory + prefixes[0] + '_' + timestamps[0] + '.h5'
    # open h5 file
    f = h5py.File(file)
    # get the dataset handle
    dset_id = f.get('blocks')

    res = dset_id.shape

    if (len(res) == 3):
        # ------------------------------------------------------
        # 2d data
        # ------------------------------------------------------
        Nb, Bs = res[0:1+1]

        # header
        fid.write('<?xml version="1.0" ?>\n')
        fid.write('<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>\n')
        fid.write('<Xdmf xmlns:xi="http://www.w3.org/2003/XInclude" Version="2.2">\n')
        fid.write('  <Domain>\n')

        fid.write('    <Grid Name="%s" GridType="Collection" CollectionType="Temporal">\n' % ('wabbit2d') )
        fid.write('\n')

        # list of times (associated to timestamps)
        fid.write('      <Time TimeType="List">\n')
        fid.write('        <DataItem Format="XML" NumberType="Float" Dimensions="%i">\n' % (len(timestamps)) )
        fid.write('          ')
        for time in times:
            fid.write(' %e ' % (time) )
        fid.write('</DataItem>\n')
        fid.write('        </Time>\n')
        fid.write('\n')

        # loop over time steps
        for i in range(len(timestamps)):
            fid.write('\n')
            fid.write('        <Grid Name="timestep %i" GridType="Collection" CollectionType="Spatial">\n' % (i))

            #use any of our files at the same timestamp to determine number of blocks
            file = directory + prefixes[0] + '_' + timestamps[i] + '.h5'
            # open h5 file
            f = h5py.File(file)
            dset_id = f.get('blocks')
            Nb = dset_id.shape[0]
            print("timestamp "+timestamps[i]+" has Nb=%i blocks" % (Nb) )

            # all blocks for this timestep
            for b  in range(Nb):
                fid.write('          <!-- ***************************************************************** -->\n')
                fid.write('          <Grid Name="block%i" GridType="Uniform">\n' % (b))
                fid.write('            <Topology TopologyType="2DCoRectMesh" NumberOfElements="%i %i"/>\n' % (Bs,Bs) )
                fid.write('            <Geometry GeometryType="ORIGIN_DXDY">\n')
                fid.write('\n')
                fid.write('              <DataItem ItemType="HyperSlab" Dimensions="2" Type="HyperSlab">\n')
                fid.write('                <DataItem Dimensions="3 2" Format="XML">\n')
                fid.write('                  %i 0\n' % (b))
                fid.write('                  1 1\n')
                fid.write('                  1 2\n')
                fid.write('                </DataItem>\n')
                fid.write('                <DataItem Name="origin" Dimensions="%i 2" Format="HDF" ItemType="Uniform">\n' % (Nb))
                # origin of blocks defined by first prefix
                fid.write('                  %s:/coords_origin\n' % (directory + prefixes[0] + '_' + timestamps[i] + '.h5') )
                fid.write('                </DataItem>\n')
                fid.write('              </DataItem>\n')
                fid.write('\n')
                fid.write('              <DataItem ItemType="HyperSlab" Dimensions="2" Type="HyperSlab">\n')
                fid.write('                <DataItem Dimensions="3 2" Format="XML">\n')
                fid.write('                  %i 0\n' % (b) )
                fid.write('                  1 1\n')
                fid.write('                  1 2\n')
                fid.write('                </DataItem>\n')
                fid.write('                <DataItem Name="spacing" Dimensions="%i 2" Format="HDF">\n' % (Nb) )
                # block spacing defined by first prefix
                fid.write('                  %s:/coords_spacing\n'% (directory + prefixes[0] + '_' + timestamps[i] + '.h5'))
                fid.write('                </DataItem>\n')
                fid.write('              </DataItem>\n')
                fid.write('            </Geometry>\n')
                fid.write('\n')
                # for each time step and each block, we have different quantities (prefixes)
                for prefix in prefixes:
                    fid.write('            <Attribute Name="%s"  AttributeType="Scalar" Center="Node">\n' % (prefix) )
                    fid.write('              <DataItem ItemType="HyperSlab" Dimensions="%i %i" Type="HyperSlab">\n' % (Bs,Bs) )
                    fid.write('                <DataItem Dimensions="3 3" Format="XML">\n')
                    fid.write('                  %i 0 0\n' % (b) )
                    fid.write('                  1 1 1\n')
                    fid.write('                  1 %i %i\n' % (Bs,Bs) )
                    fid.write('              </DataItem>\n')
                    fid.write('                <DataItem Format="HDF" Dimensions="%i %i %i" AttributeType="Scalar" Center="Node">\n' % (Nb,Bs,Bs) )
                    fid.write('                  %s:/blocks\n' % (directory + prefix + '_' + timestamps[i] + '.h5') )
                    fid.write('                </DataItem>\n')
                    fid.write('              </DataItem>\n')
                    fid.write('            </Attribute>\n')

                # although the following code appears to be correct, paraview fails to read it correctly. it reports a 3d vector
                # for the moment, we just deactivate vectors...
#                for prefix in vectors:
#                    fid.write('            <Attribute Name="%s"  AttributeType="Vector" Center="Node">\n' % (prefix) )
#                    fid.write('    <DataItem ItemType="Function" Function="JOIN($0, $1)" Dimensions="%i %i 2">\n' % (Bs, Bs) )
#                    fid.write('         <DataItem ItemType="HyperSlab" Dimensions="%i %i" Type="HyperSlab">\n' % (Bs,Bs) )
#                    fid.write('                <DataItem Dimensions="3 3" Format="XML">\n')
#                    fid.write('                  %i 0 0\n' % (b) )
#                    fid.write('                  1 1 1\n')
#                    fid.write('                  1 %i %i\n' % (Bs,Bs) )
#                    fid.write('                </DataItem>\n')
#                    fid.write('                <DataItem Format="HDF" Dimensions="%i %i %i">\n' % (Nb,Bs,Bs) )
#                    fid.write('                  %s:/blocks\n' % (directory + prefix + 'x' + '_' + timestamps[i] + '.h5') )
#                    fid.write('                </DataItem>\n')
#                    fid.write('         </DataItem>\n')
#
#                    fid.write('         <DataItem ItemType="HyperSlab" Dimensions="%i %i" Type="HyperSlab">\n' % (Bs,Bs) )
#                    fid.write('                <DataItem Dimensions="3 3" Format="XML">\n')
#                    fid.write('                  %i 0 0\n' % (b) )
#                    fid.write('                  1 1 1\n')
#                    fid.write('                  1 %i %i\n' % (Bs,Bs) )
#                    fid.write('                </DataItem>\n')
#                    fid.write('                <DataItem Format="HDF" Dimensions="%i %i %i">\n' % (Nb,Bs,Bs) )
#                    fid.write('                  %s:/blocks\n' % (directory + prefix + 'y' + '_' + timestamps[i] + '.h5') )
#                    fid.write('                </DataItem>\n')
#                    fid.write('         </DataItem>\n')
#
#                    fid.write('    </DataItem>\n')
#                    fid.write('            </Attribute>\n')

                fid.write('          </Grid>\n')
            fid.write('        </Grid> <!-- ooooooooo END OF TIME STEP oooooooo -->\n')
        fid.write('      </Grid> <!-- END OF WABBIT COLLECTION -->\n')
        fid.write('    </Domain>\n')
        fid.write('  </Xdmf>\n')
        fid.close()

    elif (len(res)==4):
        # 3d data
        Nb, Bs = res[0:1+1]
        print('File %s contains Nb=%i blocks of size %i x %i x %i)' % (file, Nb, Bs, Bs, Bs) )

        # header
        fid.write('<?xml version="1.0" ?>\n')
        fid.write('<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>\n')
        fid.write('<Xdmf xmlns:xi="http://www.w3.org/2003/XInclude" Version="2.2">\n')
        fid.write('  <Domain>\n')

        fid.write('    <Grid Name="%s" GridType="Collection" CollectionType="Temporal">\n' % ('wabbit3D') )
        fid.write('\n')

        # list of times (associated to timestamps)
        fid.write('      <Time TimeType="List">\n')
        fid.write('        <DataItem Format="XML" NumberType="Float" Dimensions="%i">\n' % (len(timestamps)) )
        fid.write('          ')
        for time in times:
            fid.write(' %e ' % (time) )
        fid.write('</DataItem>\n')
        fid.write('        </Time>\n')
        fid.write('\n')

        # loop over time steps
        for i in range(len(timestamps)):
            fid.write('\n')
            fid.write('        <Grid Name="timestep %i" GridType="Collection" CollectionType="Spatial">\n' % (i))

            #use any of our files at the same timestamp to determine number of blocks
            file = directory + prefixes[0] + '_' + timestamps[i] + '.h5'
            # open h5 file
            f = h5py.File(file)
            dset_id = f.get('blocks')
            Nb = dset_id.shape[0]
            print("timestamp "+timestamps[i]+" has Nb=%i blocks" % (Nb) )

            # all blocks for this timestep
            for b  in range(Nb):
                fid.write('          <!-- ***************************************************************** -->\n')
                fid.write('          <Grid Name="block%i" GridType="Uniform">\n' % (b))
                fid.write('            <Topology TopologyType="3DCoRectMesh" NumberOfElements="%i %i %i"/>\n' % (Bs,Bs,Bs) )
                fid.write('            <Geometry GeometryType="ORIGIN_DXDYDZ">\n')
                fid.write('\n')
                fid.write('              <DataItem ItemType="HyperSlab" Dimensions="2" Type="HyperSlab">\n')
                fid.write('                <DataItem Dimensions="3 2" Format="XML">\n')
                fid.write('                  %i 0\n' % (b))
                fid.write('                  1 1\n')
                fid.write('                  1 3\n')
                fid.write('                </DataItem>\n')
                fid.write('                <DataItem Name="origin" Dimensions="%i 3" Format="HDF" ItemType="Uniform">\n' % (Nb))
                # origin of blocks defined by first prefix
                fid.write('                  %s:/coords_origin\n' % (directory + prefixes[0] + '_' + timestamps[i] + '.h5') )
                fid.write('                </DataItem>\n')
                fid.write('              </DataItem>\n')
                fid.write('\n')
                fid.write('              <DataItem ItemType="HyperSlab" Dimensions="2" Type="HyperSlab">\n')
                fid.write('                <DataItem Dimensions="3 2" Format="XML">\n')
                fid.write('                  %i 0\n' % (b) )
                fid.write('                  1 1\n')
                fid.write('                  1 3\n')
                fid.write('                </DataItem>\n')
                fid.write('                <DataItem Name="spacing" Dimensions="%i 3" Format="HDF">\n' % (Nb) )
                # block spacing defined by first prefix
                fid.write('                  %s:/coords_spacing\n'% (directory + prefixes[0] + '_' + timestamps[i] + '.h5'))
                fid.write('                </DataItem>\n')
                fid.write('              </DataItem>\n')
                fid.write('            </Geometry>\n')
                fid.write('\n')
                # for each time step and each block, we have different quantities (prefixes)
                for prefix in prefixes:
                    fid.write('            <Attribute Name="%s"  AttributeType="Scalar" Center="Node">\n' % (prefix) )
                    fid.write('              <DataItem ItemType="HyperSlab" Dimensions="%i %i %i" Type="HyperSlab">\n' % (Bs,Bs,Bs) )
                    fid.write('                <DataItem Dimensions="3 4" Format="XML">\n')
                    fid.write('                  %i 0 0 0\n' % (b) )
                    fid.write('                  1 1 1 1 \n')
                    fid.write('                  1 %i %i %i\n' % (Bs,Bs,Bs) )
                    fid.write('              </DataItem>\n')
                    fid.write('                <DataItem Format="HDF" Dimensions="%i %i %i %i" AttributeType="Scalar" Center="Node">\n' % (Nb,Bs,Bs,Bs) )
                    fid.write('                  %s:/blocks\n' % (directory + prefix + '_' + timestamps[i] + '.h5') )
                    fid.write('                </DataItem>\n')
                    fid.write('              </DataItem>\n')
                    fid.write('            </Attribute>\n')

                # although the following code appears to be correct, paraview fails to read it correctly. it reports a 3d vector
                # for the moment, we just deactivate vectors...
#                for prefix in vectors:
#                    fid.write('            <Attribute Name="%s"  AttributeType="Vector" Center="Node">\n' % (prefix) )
#                    fid.write('    <DataItem ItemType="Function" Function="JOIN($0, $1)" Dimensions="%i %i 2">\n' % (Bs, Bs) )
#                    fid.write('         <DataItem ItemType="HyperSlab" Dimensions="%i %i" Type="HyperSlab">\n' % (Bs,Bs) )
#                    fid.write('                <DataItem Dimensions="3 3" Format="XML">\n')
#                    fid.write('                  %i 0 0\n' % (b) )
#                    fid.write('                  1 1 1\n')
#                    fid.write('                  1 %i %i\n' % (Bs,Bs) )
#                    fid.write('                </DataItem>\n')
#                    fid.write('                <DataItem Format="HDF" Dimensions="%i %i %i">\n' % (Nb,Bs,Bs) )
#                    fid.write('                  %s:/blocks\n' % (directory + prefix + 'x' + '_' + timestamps[i] + '.h5') )
#                    fid.write('                </DataItem>\n')
#                    fid.write('         </DataItem>\n')
#
#                    fid.write('         <DataItem ItemType="HyperSlab" Dimensions="%i %i" Type="HyperSlab">\n' % (Bs,Bs) )
#                    fid.write('                <DataItem Dimensions="3 3" Format="XML">\n')
#                    fid.write('                  %i 0 0\n' % (b) )
#                    fid.write('                  1 1 1\n')
#                    fid.write('                  1 %i %i\n' % (Bs,Bs) )
#                    fid.write('                </DataItem>\n')
#                    fid.write('                <DataItem Format="HDF" Dimensions="%i %i %i">\n' % (Nb,Bs,Bs) )
#                    fid.write('                  %s:/blocks\n' % (directory + prefix + 'y' + '_' + timestamps[i] + '.h5') )
#                    fid.write('                </DataItem>\n')
#                    fid.write('         </DataItem>\n')
#
#                    fid.write('    </DataItem>\n')
#                    fid.write('            </Attribute>\n')

                fid.write('          </Grid>\n')
            fid.write('        </Grid> <!-- ooooooooo END OF TIME STEP oooooooo -->\n')
        fid.write('      </Grid> <!-- END OF WABBIT COLLECTION -->\n')
        fid.write('    </Domain>\n')
        fid.write('  </Xdmf>\n')
        fid.close()




def write_xmf_file_flusi(args, outfile, times, timestamps, prefixes, scalars, vectors, directory):
    # Todo: read nx,ny,nz,lx,ly,lz,origin from files here and not in wrapper. check if every file has the same values
    # yell if not.

    print('-------------------------')
    print('- FluSI module          -')
    print('-------------------------')

    # we read the domain size and resolution from all files that we want to process
    # into lists and check if they all have the same values. If they don't, the code
    # yells but proceeds. If the resolution changes between files, the XMF file will
    # cause paraview to crash
    nnx, nny, nnz, llx, lly, llz = [], [], [], [], [], []
    ox, oy, oz = [],[],[]
    for prefix in prefixes:
        for timestamp in timestamps:
            fname = directory + prefix + '_' + timestamp + '.h5'
            # read file
            f = h5py.File(fname, 'r')

            # list all hdf5 datasets in the file - usually, we expect
            # to find only one.
            datasets = list(f.keys())

            # as there should be only one, this should be our dataset:
            dset_name = datasets[0]

            # get the dataset handle
            dset_id = f.get(dset_name)

            # from the dset handle, read the attributes
            res = dset_id.attrs.get('nxyz')
            box = dset_id.attrs.get('domain_size')
            # new: we allow for non.zero origin, if the file contains that information.
            origin = dset_id.attrs.get('origin')
            if origin is None or args.ignore_origin:
                origin = [0.0, 0.0, 0.0]

            nnx.append(res[0])
            nny.append(res[1])
            llx.append(box[0])
            lly.append(box[1])
            ox.append(origin[0])
            oy.append(origin[1])

            if len(res)==3:
                nnz.append(res[2])
                llz.append(box[2])
                oz.append(origin[2])

    print('From all files, extracting resolution and domain size...')
    print('resolution..')
    nx, ny, nz = uniquelist(nnx), uniquelist(nny), uniquelist(nnz)
    print('domain size..')
    lx, ly, lz = uniquelist(llx), uniquelist(lly), uniquelist(llz)
    print('origin..')
    ox, oy, oz = uniquelist(ox), uniquelist(oy), uniquelist(oz)
    origin = [ox,oy,oz]


    # unit spacing, if forced
    if args.unit_spacing:
        lx, ly, lz = float(nx), float(ny), float(nz)

    if len(res) == 3:
        print('Resolution is     %i x %i x %i' % (nx,ny,nz) )
        print('Domain size is    %e x %e x %e' % (lx,ly,lz) )
        print("Origin of grid is %e %e %e" % (origin[0],origin[1],origin[2]))
    else:
        print('Resolution is     %i x %i' % (nx,ny) )
        print('Domain size is    %e x %e' % (lx,ly) )
        print("Origin of grid is %e %e" % (origin[0],origin[1]))

    fid = open(outfile, 'w')

    #--------------------------------------------------------------------------
    # file header
    #--------------------------------------------------------------------------
    # NB: indices output in z,y,x order. (C vs Fortran ordering?)
    fid.write('<?xml version="1.0" ?>\n')
    fid.write('<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" [\n')
    if len(res) == 3:
        fid.write('<!ENTITY nxnynz "%i %i %i">\n' % (nz,ny,nx) )
        # write also half the resolutio to XMF file since we might use it
        # if striding is active!
        fid.write('<!ENTITY subnxnynz "%i %i %i">\n' % (nz/2, ny/2, nx/2) )
    else:
        fid.write('<!ENTITY nxnynz "%i %i">\n' % (ny,nx) )
        # write also half the resolutio to XMF file since we might use it
        # if striding is active!
        fid.write('<!ENTITY subnxnynz "%i %i">\n' % (ny/2, nx/2) )

    fid.write(']>\n')
    fid.write('<Xdmf Version="2.0">\n')
    fid.write('<Domain>\n')
    fid.write('<Grid Name="Box" GridType="Collection" CollectionType="Temporal">\n')

    for k in range(0, len(timestamps) ):
        #--------------------------------------------------------------------------
        # begin time step
        #--------------------------------------------------------------------------
        fid.write('<!-- beginning time step -->\n')
        fid.write('<Grid Name="FLUSI_cartesian_grid" GridType="Uniform">\n')
        fid.write('    <Time Value="%e" />\n' % times[k])
        fid.write('    <Topology TopologyType="%iDCoRectMesh" Dimensions="&nxnynz;" />\n' % len(res))
        fid.write(' \n')
        if len(res) == 3:
            fid.write('    <Geometry GeometryType="Origin_DxDyDz">\n')
        else:
            fid.write('    <Geometry GeometryType="Origin_DxDy">\n')
        fid.write('    <DataItem Dimensions="%i" NumberType="Float" Format="XML">\n' % len(res))
        if len(res) == 3:
            fid.write('    %e %e %e \n' % (origin[2],origin[1],origin[0]) )
        else:
            fid.write('    %e %e \n' % (origin[1], origin[0]) )
        fid.write('    </DataItem>\n')
        fid.write('    <DataItem Dimensions="%i" NumberType="Float" Format="XML">\n' % len(res))
        # NB: indices output in z,y,x order. (C vs Fortran ordering?)
        if len(res) == 3:
            fid.write('    %e %e %e\n' % (lz/nz, ly/ny, lx/nx) )
        else:
            fid.write('    %e %e\n' % (ly/ny, lx/nx) )
        fid.write('    </DataItem>\n')
        fid.write('    </Geometry>\n')

        for p in scalars:
            #++++ scalar
            # no striding, read the entire HDF5 file
            fid.write('\n')
            fid.write('    <!--Scalar-->\n')
            fid.write('    <Attribute Name="%s" AttributeType="Scalar" Center="Node">\n' % p)
            fid.write('    <DataItem Dimensions="&nxnynz;" NumberType="Float" Format="HDF">\n')
            fid.write('    %s%s_%s.h5:/%s\n' % (directory, p, timestamps[k], p) )
            fid.write('    </DataItem>\n')
            fid.write('    </Attribute>\n')

        for p in vectors:
            #+++++++++++ vector
            if len(res) == 3:
                fid.write('\n')
                fid.write('    <!--Vector-->\n')
                fid.write('    <Attribute Name="%s" AttributeType="Vector" Center="Node">\n' % p)
                fid.write('    <DataItem ItemType="Function" Function="JOIN($0, $1, $2)" Dimensions="&nxnynz; 3" NumberType="Float">\n')
                fid.write('        <DataItem Dimensions="&nxnynz;" NumberType="Float" Format="HDF">\n')
                fid.write('        %s%s_%s.h5:/%s\n' % (directory, p+'x', timestamps[k], p+'x') )
                fid.write('        </DataItem>\n')
                fid.write('\n')
                fid.write('        <DataItem Dimensions="&nxnynz;" NumberType="Float" Format="HDF">\n')
                fid.write('        %s%s_%s.h5:/%s\n' % (directory, p+'y', timestamps[k], p+'y') )
                fid.write('        </DataItem>\n')
                fid.write('\n')
                fid.write('        <DataItem Dimensions="&nxnynz;" NumberType="Float" Format="HDF">\n')
                fid.write('        %s%s_%s.h5:/%s\n' % (directory, p+'z', timestamps[k], p+'z') )
                fid.write('        </DataItem>\n')
                fid.write('    </DataItem>\n')
                fid.write('    </Attribute>\n')
            else:
                fid.write('\n')
                fid.write('    <!--Vector-->\n')
                fid.write('    <Attribute Name="%s" AttributeType="Vector" Center="Node">\n' % p)
                fid.write('    <DataItem ItemType="Function" Function="JOIN($0, $1)" Dimensions="&nxnynz; 2" NumberType="Float">\n')
                fid.write('        <DataItem Dimensions="&nxnynz;" NumberType="Float" Format="HDF">\n')
                fid.write('        %s%s_%s.h5:/%s\n' % (directory, p+'x', timestamps[k], p+'x') )
                fid.write('        </DataItem>\n')
                fid.write('\n')
                fid.write('        <DataItem Dimensions="&nxnynz;" NumberType="Float" Format="HDF">\n')
                fid.write('        %s%s_%s.h5:/%s\n' % (directory, p+'y', timestamps[k], p+'y') )
                fid.write('        </DataItem>\n')
                fid.write('    </DataItem>\n')
                fid.write('    </Attribute>\n')

        #--------------------------------------------------------------------------
        # end time step
        #--------------------------------------------------------------------------
        fid.write('</Grid>\n')
    #--------------------------------------------------------------------------
    # file footer
    #--------------------------------------------------------------------------
    fid.write('</Grid>\n')
    fid.write('</Domain>\n')
    fid.write('</Xdmf>\n')

    fid.close()
    return

def main():
    print( bcolors.OKGREEN + "**********************************************" + bcolors.ENDC )
    print( bcolors.OKGREEN + "**   hdf2xml.py                             **" + bcolors.ENDC )
    print( bcolors.OKGREEN + "**********************************************" + bcolors.ENDC )

    parser = argparse.ArgumentParser()
    parser.add_argument("-n", "--time-by-fname", help="""How shall we know at what time the file is? Sometimes, you'll end up with several
    files at the same time, which have different file names. Then you'll want to
    read the time from the filename, since paraview crashes if two files are at the
    same instant. Setting -n will force hdf2xmf.py to read from filename, eg mask_00010.h5
    will be at time 10, even if h5 attributes tell it is at t=0.1""", action="store_true")
    parser.add_argument("-1", "--one-file-per-timestep", help="""Sometimes, it is useful to generate one XMF file per
    time step (and not one global file), for example to compare two time steps. The -1 option generates these individual
    files. If -o outfile.xmf is set, then the files are named outfile_0000.xmf, outfile_0001.xmf etc.""", action="store_true")
    parser.add_argument("-o", "--outfile", help="XMF file to write to, default is ALL.xmf")
    parser.add_argument("-0", "--ignore-origin", help="force origin to 0,0,0", action="store_true")
    parser.add_argument("-u", "--unit-spacing", help="use unit spacing dx=dy=dz=1 regardless of what is specified in h5 files", action="store_true")
    parser.add_argument("-d", "--directory", help="directory of h5 files, if not ./")
    parser.add_argument("-q", "--scalars", help="""Overwrite vector recongnition. Normally, a file ux_8384.h5 is interpreted as vector,
    so we also look for uy_8384.h5 and [in 3D mode] for uz_8384.h5. -q overwrites this behavior and individually processes all prefixes as scalars.
    This option is useful if for some reason
    you have a file that ends with {x,y,z} is not a vector or if you downloaded just one component, e.g. ux_00100.h5
    """, action="store_true")
    group1 = parser.add_mutually_exclusive_group()
    group1.add_argument("-e", "--exclude-prefixes", help="Exclude these prefixes (space separated)", nargs='+')
    group1.add_argument("-i", "--include-prefixes", help="Include just these prefixes, if the files exist (space separated)", nargs='+')
    group2 = parser.add_mutually_exclusive_group()
    group2.add_argument("-t", "--include-timestamps", help="Include just use these timestamps, if the files exist (space separated)", nargs='+')
    group2.add_argument("-x", "--exclude-timestamps", help="Exclude these timestamps (space separated)", nargs='+')
    group3 = parser.add_mutually_exclusive_group()
    group3.add_argument("-p", "--skip-incomplete-timestamps", help="If some files are missing, skip the time step", action="store_true")
    group3.add_argument("-l", "--skip-incomplete-prefixes", help="If some files are missing, skip the prefix", action="store_true")
    args = parser.parse_args()

    if args.directory is None:
        directory = './'
    else:
        directory = args.directory

    if directory[-1] != "/":
        directory = directory + '/'
    print("looking for files in dir: " + bcolors.HEADER + directory + bcolors.ENDC)

    # parse the filename to write to, as in previous versions, the default value
    # is ALL.xmf
    if args.outfile == None:
        args.outfile="ALL2.xmf"
    print("XMF output will be written to: "+bcolors.HEADER+args.outfile+bcolors.ENDC)

    # How shall we know at what time the file is? Sometimes, you'll end up with several
    # files at the same time, which have different file names. Then you'll want to
    # read the time from the filename, since paraview crashes if two files are at the
    # same instant.
    if args.time_by_fname:
        print("Time will be read from: "+bcolors.HEADER + "filename" + bcolors.ENDC)
    else:
        print("Time will be read from: "+bcolors.HEADER + "dataset" + bcolors.ENDC)

    # force unit spacing, dx=dy=dz=1 (useful for micro-ct data, occasionally)
    if args.unit_spacing:
        print("We will force unit spacing! dx=dy=dz=1 regardless of what is specified in h5 files")

    if args.ignore_origin:
        print("We will force origin = 0.0, 0.0, 0.0 regardless of what is specified in h5 files")

    # will vector recognition be turned off? This option is useful if for some reason
    # you have a file that ends with x is not a vector or if you downloaded just one
    # component
    if args.scalars:
        print(bcolors.HEADER + "Vector recongnition is turned OFF! All files treated as scalars." + bcolors.ENDC)

    # it happens that you want to ignore some files with a given prefix, for example
    # if you're missing some files or want to save on memory. the --exclude option
    # lets you specify a number of space-separated prefixes for the script to ignore.
    if args.exclude_prefixes is None:
        args.exclude_prefixes = []
    print("We will exclude the following prefixes: ", end='')
    print_list(args.exclude_prefixes)

    # ...
    if args.include_prefixes is None:
        args.include_prefixes = []
    print("We will include only the following prefixes: ", end='')
    print_list(args.include_prefixes)

    # on a large dataset of files, it may be useful to ignore some time steps
    # if you're not interested in them. The --exclude-timestamps option lets you do that
    if args.exclude_timestamps is None:
        args.exclude_timestamps = []
    print("We will exclude the following timestamps: ", end='')
    print_list(args.exclude_timestamps)

    # on a large dataset of files, it may be useful to use just some time steps
    # and ignore all other.
    if args.include_timestamps is None:
        args.include_timestamps = []
    print("We will include only the following timestamps: ", end='')
    print_list(args.include_timestamps)


    # will we generate one or many XMF files?
    if args.one_file_per_timestep:
        print("XMF: One file per timestep will be written")
    else:
        print("XMF: One file with all timesteps will be generated")

    #-------------------------------------------------------------------------------
    # get the list of all h5 files in the current directory.
    #-------------------------------------------------------------------------------
    print('-------------------------------------------------------------------')
    print("Looking for files...")
    # get the list of h5 files and sort them
    filelist = sorted( glob.glob(directory + "*.h5") )
    if not filelist:
        warn('No *.h5 files found')
        return
    print("We found " + bcolors.HEADER + "%i" % len(filelist) + bcolors.ENDC + " *.h5-files in directory")

    # initialize the list of prefixes
    prefixes = []
    vectors = []
    scalars = []
    filelist_used = []
    # initialize list of times
    times = []
    # mode can be either wabbit or flusi
    mode = None

    # loop over all h5 files, add their prefix and timestamp to a list
    for file in filelist:
        # read file
        f = h5py.File(file, 'r')
        # list all hdf5 datasets in the file - usually, we expect
        # to find only one.
        datasets = list(f.keys())

        if 'blocks' in datasets and mode is None:
            # we deal with wabbit data: more than one dset per file
            mode = 'wabbit'
            print('Data identified as WABBIT data..')

        elif get_dset_name(file) in datasets and mode is None:
            # flusi data, only one dset per file expected
            mode = 'flusi'
            print('Data identified as FLUSI data..')

        if mode is None:
            # this message is also issued for FLUSI runtime backup files...
            warn('File: '+file+' seems to be neither WABBIT nor FLUSI data. Skip.')

        else:
            # prefix name
            prefix = get_dset_name(file)

            # if we find more than one dset we warn that this is unusual for flusi
            # files. for wabbit files, it is normal.
            if  (mode == 'flusi') and (len(datasets) != 1) :
                warn("we found more than one dset in the file (and thus skip it)..."+file)

            else:
                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                # judging from the dataset, do we use this file?
                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                # a priori, we'll not use this this
                used = False
                # Variant I: given list of exclude prefixes:
                if args.exclude_prefixes:
                    # the --exclude-prefixe option helps us ignore some prefixes, for example if they
                    # have incomplete data. so, if the current prefix is on that list, ignore it:
                    if prefix not in args.exclude_prefixes:
                        # we used this file:
                        used = True

                # Variant II: given list of include prefixes:
                if args.include_prefixes:
                    # the --include-prefixes option helps us focus on some prefixes, for example if they
                    # have incomplete data. so, if the current prefix is on that list, ignore it:
                    if prefix in args.include_prefixes:
                        # we used this file:
                        used = True

                # variant III: neither of both:
                if not args.exclude_prefixes and not args.include_prefixes:
                    # we used this file:
                    used = True

                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                # judging from the timestamp, do we use this file?
                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                if used:
                    # now the condition for the dataset was met, we again suppose not
                    # to use the file.
                    used = False
                    # get filename timestamp
                    timestamp = get_timestamp( file )
                    # variant I: given list of timestamps to exlude
                    if args.exclude_timestamps:
                        if timestamp not in args.exclude_timestamps:
                            # we used this file:
                            used = True

                    # variant II: given the list of timestamps
                    if args.include_timestamps:
                        if timestamp in args.include_timestamps:
                            # we used this file:
                            used = True

                    # variant III neither of both
                    if not args.exclude_timestamps and not args.include_timestamps:
                        # we used this file:
                        used = True

                if used:
                    # add to list of actually used files
                    filelist_used.append(file)
                    # store the dsetname in the list of prefixes. here, the entries are non-unique
                    # we'll remove duplicates later.
                    prefixes.append(prefix)

                    # the option --scalars forces the code to ignore the trailing x,y,z icons
                    # and treat all fields as scalars
                    # vector / scalar handling: if it ends on {x,y,z} the prefix indicates a vector
                    # otherwise, we deal with a scalar field.
                    if prefix[len(prefix)-1:len(prefix)] in ['x','y','z'] and not args.scalars:
                        # add prefix name without trailing x,y,z to list of vectors
                        vectors.append( prefix[0:len(prefix)-1] )
                    else:
                        # it's a scalar!
                        scalars.append(prefix)

    # remove duplicates
    prefixes = sorted( list(set(prefixes)) )
    vectors = sorted( list(set(vectors)) )
    scalars = sorted( list(set(scalars)) )

    #-------------------------------------------------------------------------------
    # check if vectors are complete, if not, add them to scalars (ux_00.h5 uy_00.h5 uz_00.h5)
    #-------------------------------------------------------------------------------
    for pre in vectors:
        if (pre+'x' in prefixes and pre+'y' in prefixes and pre+'z' in prefixes):
            print( pre+' is a 3D vector (x,y,z)')
        elif (pre+'x' in prefixes and pre+'y' in prefixes):
            print( pre+' is a 2D vector (x,y)')
        else:
            warn( pre+' is not a vector (its x-component is missing..)')
            vectors.remove( pre )
            if pre+'x' in prefixes:
                scalars.append(pre+'x')
            if pre+'y' in prefixes:
                scalars.append(pre+'y')
            if pre+'z' in prefixes:
                scalars.append(pre+'z')


    #-------------------------------------------------------------------------------
    # retrieve unique prefixes
    #-------------------------------------------------------------------------------
    print("We found the following prefixes: ", end='')
    print_list( prefixes )
    print("We found the following vectors: ", end='')
    print_list( vectors )
    print("We found the following scalars: ", end='')
    print_list( scalars )

    #-------------------------------------------------------------------------------
    # loop over all used files and extract timestamps
    #-------------------------------------------------------------------------------
    timestamps=[]
    for file in filelist_used:
        timestamps.append( get_timestamp( file ) )

    # retrieve unique timestamps
    timestamps = sorted( list(set(timestamps)) )
    print("We found the following timestamps: ", end='')
    print_list( timestamps )

    #-------------------------------------------------------------------------------
    # check if all files from the matrix exist
    #-------------------------------------------------------------------------------
    timestamps_to_remove, prefixes_to_remove =[], []

    for t in timestamps:
        for p in prefixes:
            # construct filename
            fname = directory + p + "_" + t + ".h5"
            if not os.path.isfile(fname):
                warn("File "+fname+ " NOT found!")

                # if desired, remove the timestamp from the list:
                if args.skip_incomplete_timestamps:
                    warn("removing timestamp "+t+ " completely!")
                    timestamps_to_remove.append(t)

                # if desired, remove the prefix from the list:
                if args.skip_incomplete_prefixes:
                    warn("removing prefix "+p+ " completely!")
                    prefixes_to_remove.append(p)
                    if not args.scalars:
                        raise ValueError("Please use --skip_incomplete_prefixes (-l) only with --scalars (-q) ")

    for t in timestamps_to_remove:
        timestamps.remove(t)

    for p in prefixes_to_remove:
        prefixes.remove(p)

    print("We found the following timestamps: ", end='')
    print_list( timestamps )

    # we have now the timestamps as an ordered list, and the times array as an ordered list
    # however, if we exclude / include some files, the lists do not match, and we select files with the
    # wrong timestamp in the xmf file.
    times=[]
    for timestamp in timestamps:
        time = None

        # if desired, we read the actual data time from the filename and not from
        # the file. it sounds silly - but it proved to be very useful, if you have two files
        # at the same time in dataset but different file name. happens not often, but happens.
        if args.time_by_fname:
            # convert the string to a float, simply.
            time = float( timestamp )
        else:
            for p in prefixes:
                fname = directory + p + "_" + timestamp + ".h5"
                # read time from file
                f = h5py.File(fname, 'r')
                # dataset name depends on program
                if mode is 'flusi':
                    dset_name = get_dset_name(fname)
                elif mode is 'wabbit':
                    dset_name = 'blocks'
                    #
                # get the dataset handle
                dset_id = f.get(dset_name)
                # from the dset handle, read the attributes
                tmp = dset_id.attrs.get('time')
                if time is None:
                    time = tmp
                else:
                    if abs(tmp-time) > 1.0e-5:
                        warn('It appears not all prefixes (with the same timestamp) are at the same time. consider using -n option. %s is at %f which is not %f' % (fname,tmp,time))

        # add time to the list of times.
        times.append( time )


    # check if the times are strictly increasing
    # PARAVIEW crashes with no clear error message if they dont
    if not strictly_increasing(times):
        print('-----------debug output----------')
        for t, tt in zip(timestamps, times):
            print("Timestamp %s time=%f" % (t,tt) )
        warn('List of times t is NOT monotonically increasing, this might cause PARAVIEW reader errors. Consider using the -n option')


    # warn if we re about to write an empty file
    if not prefixes or not timestamps:
        warn('No prefixes or timestamps..an empty file is created')

    print("The XMF file(s) refers to " + bcolors.HEADER + "%i" % (len(timestamps)*len(prefixes)) + bcolors.ENDC + " of these *.h5-files")

    if args.one_file_per_timestep:
        # extract base filename and extension
        fname, fext = os.path.splitext( args.outfile )
        # write one file per timestep
        for i in range(0, len(timestamps)):
            # construct filename
            outfile = fname + "_" + timestamps[i] + ".xmf"
            print("writing " + outfile + "....")
            if mode is 'flusi':
                write_xmf_file_flusi( args, outfile, [times[i]], [timestamps[i]], prefixes, scalars, vectors, directory)
            elif mode is 'wabbit':
                write_xmf_file_wabbit( args, outfile, [times[i]], [timestamps[i]], prefixes, scalars, vectors, directory)
    else:
        # one file for the dataset
        # write the acual xmf file with the information extracted above
        print("writing " + args.outfile + "....")
        if mode is 'flusi':
            write_xmf_file_flusi( args, args.outfile, times, timestamps, prefixes, scalars, vectors, directory)
        elif mode is 'wabbit':
            write_xmf_file_wabbit( args, args.outfile, times, timestamps, prefixes, scalars, vectors, directory)

    print("Done. Enjoy!")

# i hate python:
# LIKE, THAT IS EASY!
if __name__ == "__main__":
    main()
