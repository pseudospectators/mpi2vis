#!/bin/bash

# this is the main script to convert all *.h5 files in this directory
# to one *.xmf metafile

# note: # * this is not a "true" conversion in the sense that the
# *.xmf only tells paraview what to do with the *.h5 files. you'll
# need both!

# yes, its true! this script is completely automatical and doesn't
# need arguments

# Actually, it now takes arguments, allowing for the processing of a
# specific field.  By default, all fields are processed.  By using the
# command
#   hdf2xml.sh u
# only the velocity field (all thre components) is processed.  By
# using the command
#   hdf2xml.sh uz
# only the z-component of the velocity fields is processed.  More
# fields can be added via regexps, eg
#   hdf2xml.sh [uz,mask]
# to process only uz and the mask.  This can be extended to include
# more than two fields in the expected fashion.

# note the *.h5 files must contain the attributes "nxyz", "time",
# "domain_size". also, they must follow the file naming convention:
# mask_00010.h5 is a valid name. the dataset in the file MUST have the
# SAME name as the prefix of the file.


# Reset
Color_Off='\e[0m'       # Text Reset
# Regular Colors
Green='\e[0;32m'        # Green
Blue='\e[0;34m'         # Blue
Purple='\e[0;35m'       # Purple
Cyan='\e[0;36m'         # Cyan

echo -e $Green "**************************************" $Color_Off
echo -e $Green "**      HDF2XMF                     **" $Color_Off
echo -e $Green "**************************************" $Color_Off
#-----------------------
# Delete old files
#-----------------------
rm -f timesteps.in prefixes_vector.in prefixes_scalar.in

#-----------------------
# Find prefixes
#-----------------------
# Look through all the files whose names end with *.h5 and put them in
# the items array, as well as in the list, where file names are
# separated with colons. N is the number of items in the list.
# Exclude files containing the string "backup", which are not used
# in output files.
N=0
lastp=""
ending="h5"
if [ "$1" == "" ] ; then
    name=\*.${ending}
else
    name=${1}\*.${ending}
fi

for F in `find . -maxdepth 1 -not -name "*backup*" -name "${name}" | sort`
do
    # as is the file name with everything after
    # also remove the preceding ./ which shows up when one uses the
    # find command.
    p=$(echo ${F}  | sed 's/_[^_]*$//' | sed 's/\.\///')_
    p=${p%%_}
    if [ "$p" != "$lastp" ] ; then
	lastp=$p
	items[$N]=$p
	N=$((N+1))
    fi
done

echo -e "found prefixes: " ${items[@]} 

#-----------------------
# Indentify vectors and scalars from the prefixes
#-----------------------
# Look through all prefixed in array items[] if a prefix ends with
# "x", it's assumed to be the x-component of a vector.  the next 2
# indices then will contain "y" and "z" component. we remove them
# from the list items[] the prefix ending with "x" is then
# converted to its root (remove the trailing "x") and added to the
# list of vectors[] otherwise, we add the prefix to the list of
# scalars[]
N2=0
N3=0
for (( i=0; i<N; i++ ))
do
  # the prefix
  # echo "treating " ${items[i]}
    p=${items[i]}
    if [ "${p:${#p}-1:${#p}}" == "x" ]; then
        # is the last char an "x"? yes -> delete following entries
        unset items[i+1]
        # delete next two entrys (with y and z ending, hopefully)
        unset items[i+2]
        # the trailing "x" indicates a vector
        vectors[N2]=${p%%x} # collect entries for vector fields        
        echo ${vectors[N2]} >> ./prefixes_vector.in
        #echo "echoing "  ${vectors[N2]}
        N2=$((N2+1))
    else
    # no? it's a scalar.
        if [ "$p" != "" ]; then  # note empty values are not scalars (they are uy, uz but unset because of ux)
            scalars[N3]=${p}
            echo ${scalars[N3]} >> ./prefixes_scalar.in
            #echo "echoing "  ${scalars[N3]}
            N3=$((N3+1))           
        fi
    fi
done

if [ $N2 == 0 ]; then
    if [ $N3 == 0 ]; then
	echo "Error: no input data found. Exiting."
	exit
    fi
fi

# Print summary
echo "Number of vectors: "$N2
echo "Number of scalars: "$N3
echo -e "found scalars : " ${Cyan} ${scalars[@]} ${Color_Off}
echo -e "found vectors : " ${Cyan} ${vectors[@]} ${Color_Off}


# Look for time steps
if [ $N3 != 0 ]; then
    # echo "Look for time with scalars."
    i=0
    for F in `ls ${scalars[0]}*.${ending}`
    do
	time=${F%%.${ending}}
	time=${time##${scalars[0]}}
	time=${time##_}
	
	all_times[i]=${time}
	i=$((i+1))
	
	echo $time >> ./timesteps.in
    done
else 
    if [ $N2 != 0 ]; then
	# echo "Look for time with vectors."
	i=0
	for F in `ls ${vectors[0]}x*.${ending}`
	do
	    time=${F%%.${ending}}
	    time=${time##${vectors[0]}}
	    time=${time##x_}
	    
	    all_times[i]=${time}
	    i=$((i+1))
	    
	    echo $time >> ./timesteps.in
    done
    fi
fi

echo -e "found times :   " ${Cyan} ${all_times[@]} ${Color_Off}

# Create All.xmf using the FORTRAN converter
convert_hdf2xmf

# Remove temporary files.
rm -f timesteps.in prefixes_vector.in prefixes_scalar.in
