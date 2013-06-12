#!/bin/bash
#-------------------------------------------------------------------------------
# this is the main script to convert all *.h5 files in this directory to one *.xmf metafile
# note:
# * this is not a "true" conversion in the sense that the *.xmf only tells paraview what
#   to do with the *.h5 files. you'll need both!
# * yes, its true! this script is completely automatical and doesn't need arguments
# * note the *.h5 files must contain the attributes "nxyz", "time", "domain_size". also, they must
#   follow the file naming convention: mask_00010.h5 is a valid name. the dataset in the file MUST
#   have the SAME name as the prefix of the file.
#-------------------------------------------------------------------------------


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
echo -e $Green "**      T. Engels                   **" $Color_Off
echo -e $Green "**      Aix-Marseille-Université    **" $Color_Off
echo -e $Green "**************************************" $Color_Off


#-----------------------
# delete old files
#-----------------------
rm timesteps.in prefixes_vector.in prefixes_scalar.in


#-----------------------
# find prefixes
#-----------------------
# look through all the files whose names end with *.h5 and put
# them in the items array, as well as in the list, where file names
# are separated with colons. N is the number of items in the list.
N=0
lastp=""
ending="h5"
for F in `ls *.${ending}`
do
    # as is the file name with everything after
    p=$(echo ${F}  | sed 's/_[^_]*$//')_
    p=${p%%_}
    if [ "$p" != "$lastp" ] ; then
	lastp=$p
	items[$N]=$p
	N=$((N+1))
    fi
done

echo -e "found prefixes: " ${Blue} ${items[@]} ${Color_Off}



#-----------------------
# indentify vectors and scalars from the prefixes
#-----------------------
# look through all prefixed in array items[] if a prefix ends with
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
    p=${items[i]}
    if [ "${p:${#p}-1:${#p}}" == "x" ]; then
        # is the last char an "x"? yes -> delete following entries
        unset items[i+1]
        # delete next two entrys (with y and z ending, hopefully)
        unset items[i+2]
        # the trailing "x" indicates a vector
        vectors[N2]=${p%%x} # collect entries for vector fields        
        echo ${vectors[n2]} >> prefixes_vector.in
        N2=$((N2+1))
    else
    # no? it's a scalar.
        if [ "$p" != "" ]; then  # note empty values are not scalars (they are uy, uz but unset because of ux)
            scalars[N3]=${p}
            echo ${scalars[n2]} >> prefixes_scalar.in
            N3=$((N3+1))           
        fi
    fi
done

# print summary
echo -e "found scalars : " ${Cyan} ${scalars[@]} ${Color_Off}
echo -e "found vectors : " ${Cyan} ${vectors[@]} ${Color_Off}


i=0
for F in `ls ${scalars[0]}*.${ending}`
do
    time=${F%%.${ending}}
    time=${time##${scalars[0]}}
    time=${time##_}
    
    all_times[i]=${time}
    i=$((i+1))
    
    echo $time >> timesteps.in
done

echo -e "found times :   " ${Cyan} ${all_times[@]} ${Color_Off}

echo -e ${Purple} "any key to continue" ${Color_Off}
read dummy




convert_hdf2xmf





