#!/bin/bash

if [ "$1" = "-hht" ]; then
	#=================================================================
	# Hilbert Huang Transform
	#=================================================================
	read -p "Input file name: " MAPI
	read -p "Number of Intrinsic Mode Functions: " IMF
	read -p "Number of iterations for each IMF: " NIT
	read -p "Stiffness parameter: " STF
    read -p "Tension parameter: " TNS
    read -p "FWHM parameter: " GSP
	time ./hht $MAPI $IMF $NIT $STF $TNS $GSP
	#=================================================================
elif [ "$1" = "-extrema" ]; then
	#=================================================================
	# Get locsl extrema
	#=================================================================
	# Parameters to obtain local extrema maps
	read -p "Input file name: " MAPI
	read -p "Stiffness parameter: " STF
    read -p "Tension parameter: " TNS
    read -p "FWHM parameter: " GSP
	time ./hht -lext $MAPI $STF $TNS $GSP
	#=================================================================
elif [ "$1" = "-images" ]; then
	#=================================================================
	# Images
	#=================================================================
	# Parameters for MAP2GIF
	M2G=/Users/guillermo/Documents/Healpix_3.82/bin/map2gif
	RES=1900
	BAR=.true.
	read -p "Kind of map (CMB=1 or DUST=2): " MAPK
	# Choose image parameters
	if [ $MAPK = 1 ]; then
		#CMB
		LOG=.false.
		MIN=-5E-4
		MAX=5E-4
	elif [ $MAPK = 2 ]; then
		#DUST
		LOG=.true.
		MIN=-7
		MAX=-1
	else
		echo "Cannot parse argument."
	fi
	read -p "Apply default bounds for colorbar (yes=1 or no=2): " MNMX
	# Generate images
	if [ $MNMX = 1 ]; then
		for MAP in {imf*.fits,res*.fits,emd.fits,loc_*.fits,int_*.fits}
		do
			NAME="${MAP%.*}"
			$M2G -inp $NAME.fits -out $NAME.gif -xsz $RES -bar $BAR -log $LOG -min $MIN -max $MAX
		done
	elif [ $MNMX = 2 ]; then
		for MAP in {imf*.fits,res*.fits,emd.fits,loc_*.fits,int_*.fits}
		do
			NAME="${MAP%.*}"
			$M2G -inp $NAME.fits -out $NAME.gif -xsz $RES -bar $BAR -log $LOG
		done
	else
		echo "Cannot parse argument."
	fi
	#=================================================================
elif [ "$1" = "-clean" ]; then
	#=================================================================
	# Remove files
	#=================================================================
	read -p "Files to remove (FITS=1, GIF=2 or all=3): " RMV
	if [ $RMV = 1 ]; then
		for MAP in {imf*.fits,res*.fits,emd.fits,loc_*.fits,int_*.fits}
		do
			rm -fv $MAP
		done
	elif [ $RMV = 2 ]; then
		for MAP in {imf*.fits,res*.fits,emd.fits,loc_*.fits,int_*.fits}
		do
			NAME="${MAP%.*}"
			rm -fv $NAME.gif
		done
	elif [ $RMV = 3 ]; then
		for MAP in {imf*.fits,res*.fits,emd.fits,loc_*.fits,int_*.fits}
		do
			rm -fv $MAP
			NAME="${MAP%.*}"
			rm -fv $NAME.gif
		done
	else
		echo "Cannot parse argument."
	fi
	#=================================================================
else
	echo "Cannot parse 1st argument."
fi

