#!/bin/bash

# this is the main script to convert all *.h5 files in this directory
# to one *.xmf metafile

# note: # * this is not a "true" conversion in the sense that the
# *.xmf only tells paraview what to do with the *.h5 files. you'll
# need both!

# yes, its true! this script is completely automatical and doesn't
# need arguments

# The script hdf2xml.sh can be run without arguments as usual.
# By default, all fields are processed.  By using the command
#
# hdf2xml.sh -i u
#
# only the velocity field (all three components) is processed.  By
# using the command
#
# hdf2xml.sh -i uz
#
# only the z-component of the velocity fields is processed.
# using the command
#
# hdf2xml.sh -i uz -e mask
#
# only the z-component of the velocity fields is processed, but the mask
# isn't
#
# hdf2xml.sh -e mask -s
#
# all fields are processed, except the mask.
# Also, load only every 2nd grid point ("-s").
#
# ./hdf2xml.sh -i [u,mask] -e us
#
# processes u, full vector, but not us, also mask.
#
# hdf2xml.sh -p "mask,iy,usz"
#
# uses only the list of comma separated prefixes. do not forget the quotes!


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

echo -e $Green "*************************************************" $Color_Off
echo -e $Green "**      HDF2XMF                                **" $Color_Off
echo -e $Green "** ./hdf2xml.sh -i [INCLUDE] -e [EXCLUDE] -s   **" $Color_Off
echo -e $Green "**  -s Striding. Use only every 2nd grid point **" $Color_Off
echo -e $Green "**  -i prefixes to include                     **" $Color_Off
echo -e $Green "**  -e prefixes to exclude                     **" $Color_Off
echo -e $Green "**  -p specifiy prefixed manually              **" $Color_Off
echo -e $Green "**  -o specifiy outfile                        **" $Color_Off
echo -e $Green "*************************************************" $Color_Off
echo -e $Green "** Read full files u, but not mask:            **" $Color_Off
echo -e $Green "** ./hdf2xml.sh -i u -e mask                   **" $Color_Off
echo -e $Green "*************************************************" $Color_Off
echo -e $Green "** Read full files u,mask, but not us:         **" $Color_Off
echo -e $Green "** ./hdf2xml.sh -i [u,mask] -e us              **" $Color_Off
echo -e $Green "*************************************************" $Color_Off
echo -e $Green "** Read full files u,mask, but not us:         **" $Color_Off
echo -e $Green "** ./hdf2xml.sh -p \"ux,uy,uz,mask\"             **" $Color_Off
echo -e $Green "** do not forget to set the quotes !!          **" $Color_Off
echo -e $Green "*************************************************" $Color_Off

prefixes_in=""
outfile="ALL.xmf"

# parse options
while getopts ':se:i:p:o:' OPTION ; do
  case "$OPTION" in
    s)   echo -e ${Purple} "Use striding!" ${Color_Off} ; stride="y";;
    i)   echo -e "Include the following files" ${Purple} ${OPTARG} ${Color_Off} ; include=${OPTARG};;
    e)   echo -e "Exclude the following files" ${Purple} ${OPTARG} ${Color_Off} ; exclude=${OPTARG};;
    p)   echo -e "Use only prefixes" ${Purple} ${OPTARG} ${Color_Off} ; prefixes_in=${OPTARG};;
    o)   echo -e "write to outfile" ${Purple} ${OPTARG} ${Color_Off} ; outfile=${OPTARG};;
    *)   echo "Unknown parameter" ; exit 1 ;;
  esac
done

#-----------------------
# Delete old files
#-----------------------
rm -f timesteps.in prefixes_vector.in prefixes_scalar.in STRIDE.in

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


# the call to the script may predictate what prefixes to use
if [ "$prefixes_in" != "" ] ; then
  # the list of prefixes with the -p option is a string
  # with comma separated values. now we convert this string into an array:
  IFS=',' read -r -a items <<< "$prefixes_in"
  N=${#items[@]}
else
  # no list of prefixes is given, so we look for them using the find command.
  # it is possibly modified by the -e and -i options
  # explicitly include these files
  if [ "$include" == "" ] ; then
    names_include=\*.${ending}
  else
    names_include=$include\*.${ending}
  fi

  # explicitly EXCLUDE these files:
  if [ "$exclude" == "" ] ; then
    names_exclude=""
  else
    names_exclude=$exclude\*.${ending}
  fi


  for F in `find . -maxdepth 1 -not -name "*backup*" -name "${names_include}" -not -name "${names_exclude}" | sort`
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
fi

echo -e "found prefixes: " ${Cyan} ${items[@]} ${Color_Off}
all_prefixes=${items[@]}

if [ "$stride" == "y" ] ; then
  # the fortran program will check if the file STRIDE.in
  # exists and if it does, then we use striding
  echo "STRIDE ON"
  touch STRIDE.in
fi

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
echo -e "Number of vectors: "${Cyan} $N2 ${Color_Off}
echo -e "Number of scalars: "${Cyan} $N3 ${Color_Off}
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


# check if some files do not exist and if so, yell at user.
for prefix in $all_prefixes
do
  for tt in ${all_times[@]}
  do
    if [ ! -f ${prefix}"_"${tt}".h5" ]; then
      echo -e ${Cyan} file not found: ${prefix}"_"${tt}".h5"  ${Color_Off}
    fi
  done
done



echo "Do you want to create one ALL.xmf file with all time steps or one xmf-file for each time step?"
echo "[return] for ALL.xmf, (i) for individual files"
read all

if [ "$all" == "i" ]; then
  # by choice, we create one XMF file PER time step. This is handy for very large data sets.
  for time in ${all_times[@]}
  do
    rm -f timesteps.in
    echo $time >> ./timesteps.in
    # Create All.xmf using the FORTRAN converter
    convert_hdf2xmf
    mv ALL.xmf $time.xmf
  done
else
  # Create All.xmf using the FORTRAN converter
  convert_hdf2xmf
fi

# rename file to what is desired (mostly ALL.xmf)
# not ALL.xmf is erased in any case
mv ALL.xmf $outfile

# Remove temporary files.
rm -f timesteps.in prefixes_vector.in prefixes_scalar.in STRIDE.in
