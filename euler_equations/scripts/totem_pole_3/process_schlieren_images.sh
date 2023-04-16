#!/bin/bash 
set -e
set -x

count=0
name="totem_pole_3"
out_dir=./../../artwork_output/${name}
mkdir -p $out_dir
mkdir -p ${out_dir}/pdfs
mkdir -p ${out_dir}/pngs

# Iterate over all schlieren .dat files
for dat_file in ./../../euler_simulation_code/output/*schlieren*.dat; do

    out_file=${name}_frame_${count}

    echo "-----"
    echo "Processing $dat_file into $out_file:"

    
    # Run the Gnuplot script with the current .dat file as an input
    gnuplot -e "output_file='${out_dir}/${out_file}.eps'; data_file='${dat_file}'" plot_schlieren.gp
    echo "Created eps using gnuplot"
    
    # Convert the EPS plot to a PDF plot
    epstopdf "${out_dir}/${out_file}.eps" --output="${out_dir}/pdfs/${out_file}.pdf"
    echo "Converted to pdf"

    convert -density 200 "${out_dir}/pdfs/${out_file}.pdf" -quality 100 "${out_dir}/pngs/${out_file}.png"
    echo "Converted to png"
    
    # Delete the EPS plot
    echo "Deleting eps file ${out_file}.eps" && rm ${out_dir}/${out_file}.eps
    echo "-----"

    let count=count+1


    echo "Count=$count"

done

# Finally create a gif from the pns
mkdir ${out_dir}/gif
convert -delay 75 -loop 0 $(ls ${out_dir}/pngs/*_*.png | sort -V) ${out_dir}/gif/${out_file}.gif

echo "Deleting all state output files..." && rm ./../../euler_simulation_code/output/*state*.dat
echo "Deleting all schlieren output files..." && rm ./../../euler_simulation_code/output/*schlieren*.dat


