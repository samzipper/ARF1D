#!/bin/bash

for i in {1..100}
		do
		# define the prefix
		num=`printf %04d $i`

		# copy the appropriate geotop input file, which is created with R
		cp $(pwd)/SoilPropertyCalibration/${num}_SoilNRCS001.txt $(pwd)/soil/SoilNRCS0001.txt

		# run geotop
		docker run --rm -v $(pwd):/work omslab/geotop

		# copy output file
		cp $(pwd)/output-tabs/soiltemp0001.txt $(pwd)/SoilPropertyCalibration/output_${num}_soiltemp0001.txt
		cp $(pwd)/output-tabs/thetaliq0001.txt $(pwd)/SoilPropertyCalibration/output_${num}_thetaliq0001.txt

		# status update
		echo "$num" complete
	done
