# Check status of running "SLURMed" URANOS jobs
# Prints the last line of all output files in a folder.

# 1. Specify output folder here.
#    Should contain *.out files created by the SLURM jobs
for file in /work/%u/URANOS/log/URANOS-*
do
   lastline=$(sed -n '/^[[:space:]]*$/!h;$!d;g;s/^/   100 % : /p' "$file")
   echo "$lastline $file"
done

# 2. Call
#    Change chmod to 0777 and call with ./check_uranos.sh