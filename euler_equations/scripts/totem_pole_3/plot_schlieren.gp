set terminal postscript eps enhanced
set output output_file
set size 12,12
set size square
unset xlabel
unset ylabel
unset title
unset colorbox
unset border
unset tics

plot data_file u 1:2:(-exp(-20 * $3)) w image, data_file u (- $1 + 1.0/6144.0):2:(-exp(-20 * $3)) w image
