#!/bin/python
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

def print_list( l ):
    for p in l:
        print(bcolors.HEADER + p + " " + bcolors.ENDC, end='')
    # print just one newline
    print(' ')

def warn( msg ):
    print( bcolors.FAIL + "WARNING! " + bcolors.ENDC + msg)


def write_xmf_file(outfile,nx,ny,nz,lx,ly,lz, times, timestamps, prefixes, scalars, vectors, dims, directory):
    fid = open(outfile, 'w')

    #--------------------------------------------------------------------------
    # file header
    #--------------------------------------------------------------------------
    # NB: indices output in z,y,x order. (C vs Fortran ordering?)
    fid.write('<?xml version="1.0" ?>\n')
    fid.write('<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" [\n')
    if dims == 3:
        fid.write('<!ENTITY nxnynz "%i %i %i">\n' % (nz,ny,nx) )
        # write also half the resolutio to XMF file since we might use it
        # if striding is active!
        fid.write('<!ENTITY subnxnynz "%i %i %i">\n' % (nz/2, ny/2, nx/2) )
    else:
        fid.write('<!ENTITY nxnynz "%i %i">\n' % (nz,ny) )
        # write also half the resolutio to XMF file since we might use it
        # if striding is active!
        fid.write('<!ENTITY subnxnynz "%i %i">\n' % (nz/2, ny/2) )

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
        fid.write('    <Topology TopologyType="%iDCoRectMesh" Dimensions="&nxnynz;" />\n' % dims)
        fid.write(' \n')
        if dims == 3:
            fid.write('    <Geometry GeometryType="Origin_DxDyDz">\n')
        else:
            fid.write('    <Geometry GeometryType="Origin_DxDy">\n')
        fid.write('    <DataItem Dimensions="%i" NumberType="Float" Format="XML">\n' % dims)
        if dims == 3:
            fid.write('    0 0 0\n')
        else:
            fid.write('    0 0\n')
        fid.write('    </DataItem>\n')
        fid.write('    <DataItem Dimensions="%i" NumberType="Float" Format="XML">\n' % dims)
        # NB: indices output in z,y,x order. (C vs Fortran ordering?)
        if dims == 3:
            fid.write('    %e %e %e\n' % (lz/nz, ly/ny, lx/nx) )
        else:
            fid.write('    %e %e\n' % (lz/nz, ly/ny) )
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
            if dims == 3:
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
                fid.write('        %s%s_%s.h5:/%s\n' % (directory, p+'y', timestamps[k], p+'x') )
                fid.write('        </DataItem>\n')
                fid.write('\n')
                fid.write('        <DataItem Dimensions="&nxnynz;" NumberType="Float" Format="HDF">\n')
                fid.write('        %s%s_%s.h5:/%s\n' % (directory, p+'z', timestamps[k], p+'y') )
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
    parser.add_argument("-2", "--two-dim", help="Assume 2D data", action="store_true")
    parser.add_argument("-1", "--one-file-per-timestep", help="""Sometimes, it is useful to generate one XMF file per
    time step (and not one global file), for example to compare two time steps. The -1 option generates these individual
    files. If -o outfile.xmf is set, then the files are named outfile_0000.xmf, outfile_0001.xmf etc.""", action="store_true")
    parser.add_argument("-o", "--outfile", help="XMF file to write to, default is ALL.xmf")
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
    args = parser.parse_args()

    if args.directory is None:
        directory = './'
    else:
        directory = args.directory
    print("looking for files in dir: " + bcolors.HEADER + directory + bcolors.ENDC)

    # check if we deal with2d or 3d data
    print(args)
    if args.two_dim:
        dims = 2
    else:
        dims = 3
    print("HDF5 files will be treated as: "+bcolors.HEADER+("%iD data" % dims)+bcolors.ENDC)

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
    print("Looking for files...", end='')
    # get the list of h5 files and sort them
    filelist = sorted( glob.glob(directory + "*.h5") )
    if not filelist:
        warn('No files found')
        return

    # initialize the list of prefixes
    prefixes = []
    vectors = []
    scalars = []
    filelist_used = []
    # initliazie list of times
    times = []
    nx = None
    ny = None
    nz = None
    nfiles = 0

    # loop over all h5 files, add their prefix and timestamp to a list
    for file in filelist:
        # read file
        f = h5py.File(file, 'r')
        # list all hdf5 datasets in the file - usually, we expect
        # to find only one.
        datasets = f.keys()

        # if we find more than one dset we warn that this is unusual
        if (len(datasets) != 1):
            warn("we found more than one dset in the file (and thus skip it)..."+file)
        else:
            # as there should be only one, this should be our dataset:
            dset_name = datasets[0]

            # a priori, we'll not use this this
            used = False

            # Variant I: given list of exclude prefixes:
            if args.exclude_prefixes:
                # the --exclude-prefixe option helps us ignore some prefixes, for example if they
                # have incomplete data. so, if the current dset_name is on that list, ignore it:
                if dset_name not in args.exclude_prefixes:
                    # we used this file:
                    used = True

            # Variant II: given list of include prefixes:
            if args.include_prefixes:
                # the --include-prefixes option helps us focus on some prefixes, for example if they
                # have incomplete data. so, if the current dset_name is on that list, ignore it:
                if dset_name in args.include_prefixes:
                    # we used this file:
                    used = True

            # variant II: neither of both:
            if not args.exclude_prefixes and not args.include_prefixes:
                # we used this file:
                used = True

            if used:
                # add to list of actually used files
                filelist_used.append(file)
                # count used files
                nfiles = nfiles + 1
                # store the dsetname in the list of prefixes. here, the entries are non-unique
                # we'll remove duplicates later.
                prefixes.append(dset_name)

                # get the dataset handle
                dset_id = f.get(dset_name)

                # from the dset handle, read the attributes
                time = dset_id.attrs.get('time')
                times.append(time)
                res = dset_id.attrs.get('nxyz')
                box = dset_id.attrs.get('domain_size')

                # copy current values (these should not change between files. TODO: warn if they do)
                if nx is not None:
                    if nx != res[0]:
                        warn(' The nx-resolution seems to have changed')
                if ny is not None:
                    if ny != res[1]:
                        warn(' The ny-resolution seems to have changed')
                if nz is not None:
                    if nz != res[2]:
                        warn(' The nz-resolution seems to have changed')

                nx = res[0]
                ny = res[1]
                nz = res[2]
                lx = box[0]
                ly = box[1]
                lz = box[2]

                # warn if we might wrongly treat 2D data:
#                if nx == 1 and dims==3:
#                    warn('nz==1, so you might consider setting the -2 option')

                # the option --scalars forces the code to ignore the trailing x,y,z icons
                # and treat all fields as scalars
                # vector / scalar handling: if it ends on {x,y,z} the dset_name indicates a vector
                # otherwise, we deal with a scalar field.
                if dset_name[len(dset_name)-1:len(dset_name)] in ['x','y','z'] and not args.scalars:
                    vectors.append( dset_name[0:len(dset_name)-1] )
                else:
                    scalars.append(dset_name)


    # what is the resolution?
    print("Resolution is %i %i %i" % (nx,ny,nz))

    prefixes = sorted( list(set(prefixes)) )
    vectors = sorted( list(set(vectors)) )
    scalars = sorted( list(set(scalars)) )
    #-------------------------------------------------------------------------------
    # check if vectors are complete, if not, add them to scalars (ux_00.h5 uy_00.h5 uz_00.h5)
    #-------------------------------------------------------------------------------
    for pre in vectors:
        if pre+'x' not in prefixes or pre+'y' not in prefixes or pre+'z' not in prefixes:
            warn(pre+' is not a vector')
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
    # loop over all files and extract timestamps
    #-------------------------------------------------------------------------------
    import re
    timestamps=[]
    for file in filelist_used:
        # extract everything between "_" and "." so mask_00000.h5 gives 00000
        # note escape character "\"
        m = re.search('\_(.+?)\.', os.path.basename(file) )
        if m:
            found = m.group(1)
            # variant I: given list of timestamps to exlude
            if args.exclude_timestamps:
                if found not in args.exclude_timestamps:
                    # append the extracted time stampt to the list of timestamps
                    timestamps.append(found)

            # variant II: given the list of timestamps
            if args.include_timestamps:
                if found in args.include_timestamps:
                    # append the extracted time stampt to the list of timestamps
                    timestamps.append(found)

            # variant II neither of both
            if not args.exclude_timestamps and not args.include_timestamps:
                if found:
                    timestamps.append(found)
        else:
            print("An error occured: we couldn't extract the timestamp")

    # retrieve unique timestamps
    timestamps = sorted( list(set(timestamps)) )
    print("We found the following timestamps: ", end='')
    print_list( timestamps )


    # if desired, we read the actual data time from the filename and not from
    # the file. it sounds silly - but it proved to be very useful, if you have two files
    # at the same time in dataset but different file name. happens not often, but happens.
    if args.time_by_fname:
        for i in range(0, len(timestamps) ):
            # convert the string to a float, simply.
            times[i] = float( timestamps[i] )

    # as a last step before XMF generation, check if all files from the matrix exist
    for t in timestamps:
        for p in prefixes:
            # construct filename
            fname = directory + p + "_" + t + ".h5"
            if not os.path.isfile(fname):
                warn("File "+fname+ " NOT found!")

    # we have now the timestamps as an ordered list, and the times array as an ordered list
    # however, if we exlcude / include some files, the lists do not match, and we select files with the
    # wrong timestamp in the xmf file.
    # So now, we just take one prefix, loop over all used timestamps, and read the time from
    # that file. then we're sure both match.
    p = prefixes[0]
    times=[]
    for t in timestamps:
        # read file
        f = h5py.File(directory + p + "_" + t + ".h5", 'r')
        # list all hdf5 datasets in the file - usually, we expect
        # to find only one.
        datasets = f.keys()
        # as there should be only one, this should be our dataset:
        dset_name = datasets[0]
        # get the dataset handle
        dset_id = f.get(dset_name)
        # from the dset handle, read the attributes
        times.append( dset_id.attrs.get('time') )

    # warn if we re about to write an empty file
    if not prefixes or not timestamps:
        warn('No prefixes or timestamps..an empty file is created')

    print("We found " + bcolors.HEADER + "%i" % len(filelist) + bcolors.ENDC + " *.h5-files in directory")
    print("The XMF file(s) refers to " + bcolors.HEADER + "%i" % (len(timestamps)*len(prefixes)) + bcolors.ENDC + " of these *.h5-files")

    if args.one_file_per_timestep:
        # extract base filename and extension
        fname, fext = os.path.splitext( args.outfile )
        # write one file per timestep
        for i in range(0, len(timestamps)):
            # construct filename
            outfile = fname + "_" + timestamps[i] + ".xmf"
            print("writing " + outfile + "....")
            write_xmf_file( outfile,nx,ny,nz,lx,ly,lz, [times[i]], [timestamps[i]], prefixes, scalars, vectors, dims, directory )
    else:
        # one file for the dataset
        # write the acual xmf file with the information extracted above
        print("writing " + args.outfile + "....")
        write_xmf_file( args.outfile,nx,ny,nz,lx,ly,lz, times, timestamps, prefixes, scalars, vectors, dims, directory )

    print("Done. Enjoy!")

# i hate python:
# LIKE, THAT IS EASY!
if __name__ == "__main__":
    main()
