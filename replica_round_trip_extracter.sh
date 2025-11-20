#!/bin/bash

# whoami
echo "# This script reads a .log file from a GROMACS simulation with replica exchange and returns the number of times and time periods or replica 0 returning back into replica 0"

# check if we have inputs
if [ $# -lt 1 ]; then
    echo "Usage: $0 <file.log> [-histo #nbins]"
    echo "Option no histogram: $0 <file.log>"
    echo "Option histogram with X bins: $0 <file.log> -histo X"
    exit 1
fi

# Check is files exists
if [ -z "$1" ]; then
    echo "Error: No input file provided."
    exit 1
fi

# Check for optional flag for histogram generation
if [[ "$2" == "-histo" ]]; then
    histo_flag=1
elif [[ "$2" != "-histo" && ${#2} -gt 0 ]]; then
    echo "Error: only flag supported is -histo for histogram generation, provided flag "${2}" doesn't exist."
    exit 1
else
    histo_flag=0
fi

# If histogram is selected, check that a feasible number of bins is requested
if [[ ${histo_flag} -eq 1 ]]; then
    if [[ "$3" -gt 5 ]] && [[ ${#3} -gt 0 ]]; then
        bins=$3
    else
        if [[ ${#3} -eq 0 ]]; then bins=0;else bins=$3;fi
	echo "Error: you are requesting "${bins}" bins, but it should be at least 5 (and positive integer...)."
        exit 1
    fi
fi

# Prepare output filename
output="replica_time_${1%.*}.dat"

# Gather some data from the log file and prepare for replica analysis
echo "# Input file: "$1 | tee ${output}
day=$(date "+%d-%m-%Y")
time=$(date "+%Hh%Mm%Ss")
echo "# Run on "${day}" at "${time} >> ${output}

nrepl=$(grep -m 1 "Repl  There are" $1 | awk '{print $4}')
echo "# Number of replicas: "${nrepl} | tee -a ${output}

if [ ${nrepl} -lt 3 ]; then
 echo "# This number of replicas is not supported."
 exit 1
fi

temp=$(grep -m 1 "-replex" $1)
nsteps=$(echo ${temp} | awk '{for(i=1;i<=NF;i++) if ($i=="-replex") print $(i+1)}')
dt=$(grep -m 1 " dt " $1 | awk '{print $3}')
deltat=$(echo "${nsteps}*${dt}" | bc)
echo "# Simulation timestep: "${dt}" ps" | tee -a ${output}
echo "# Replica exchange was attempted every "${nsteps}" steps ("${deltat}" ps)" | tee -a ${output}

t0=$(grep -m 1 "Replica exchange at step " $1 | awk '{print $7}')
t0=$(echo "${t0}-${deltat}" | bc)
echo "# Simulation started at time "${t0}" ps" | tee -a ${output}

# Check replica exchanges
awk -v nrepl=$((nrepl-1)) -v t0=${t0} -v deltat=${deltat} -v output=${output} '
BEGIN { iamin=0;end=0;leave_time=t0 }
# Returns (first) field in a line where val == i
function find_replica(val,i)
{
    for (i=1;i<=NF;i++) if ($i == val) return i
}

{
    # Consider only replica exchange lines of log file
    if ($1" "$2 == "Repl ex")
    {
	# Update time
        frames++
	time = t0 + frames*deltat

	# gets field of where the replica 0 is
        pos = find_replica(iamin)
        
        # forward move
        if ($(pos+1) == "x")
	{
	    iamin++
	    if (iamin == nrepl && end != 1) 
            {
            end=1
            }
            if (iamin == 1) 
            {
            leave_time = time - deltat
            }
	}
        
	# backward move
        if ($(pos-1) == "x")
        {
            iamin--
	    if (iamin == 0 && end == 1)
            {
	    return_time = time - leave_time
	    round_trips++
	    # Store leave, return, and trip times on completed round trip
	    leave_times[round_trips] = leave_time/1000
            arrive_times[round_trips] = time/1000
	    trip_times[round_trips] = return_time/1000
            end=0
            }
        }
    }
}
END{
    # Summarize output of analysis
    sim_length = time + deltat
    print "# Total simulation time: " sim_length/1000 " ns"
    print "# Total simulation time: " sim_length/1000 " ns" >> output
    print "# Total round trips of replica 0: " round_trips
    print "# Total round trips of replica 0: " round_trips >> output

    # Calculate average and std dev
    ave_trip_time = 0
    for (i = 1; i <= round_trips; i++) {
        ave_trip_time += trip_times[i]
    }
    ave_trip_time /= round_trips
    std_dev_trip_time = 0
    for (i = 1; i <= round_trips; i++) {
        std_dev_trip_time += (ave_trip_time - trip_times[i])^2
    }
    std_dev_trip_time = sqrt(std_dev_trip_time/(round_trips-1))

    print "# Average round trip time " ave_trip_time" ns"
    print "# Average round trip time " ave_trip_time" ns" >> output
    print "# Std dev round trip time " std_dev_trip_time" ns"
    print "# Std dev round trip time " std_dev_trip_time" ns" >> output
    print "# Data collected in "output

    print "# Round trip times are in ns" >> output
    print "# Round_trip_number leave_time arrive_time tot_time" >> output
    for (i=1; i<=round_trips; i++) {
        print i" "leave_times[i]" "arrive_times[i]" "trip_times[i] >> output
    }
}
' $1

# Do you want a fancy histogram?
if [ ${histo_flag} -eq 1 ]; then
    # Read the trip times from the produced file and get min, max, and bin width
    mapfile -t trip_times < <(awk '{if ($1 != "#") {print $4}}' "${output}")
    #bins=20
    min=$(printf '%s\n' "${trip_times[@]}" | sort -n | head -1)
    max=$(printf '%s\n' "${trip_times[@]}" | sort -n | tail -1)
    bin_width=$(echo "(${max}-${min}) / ${bins}" | bc -l)

    # Initialize histogram
    declare -a hist
    for ((i=0; i<bins; i++)); do hist[i]=0; done

    # Fill histogram
    for time in "${trip_times[@]}"
    do
        bin=$(echo "scale=0; (${time} - ${min}) / ${bin_width}" | bc)    
        if [ ${bin} -ge ${bins} ]; then bin=$((bins - 1)); fi    
        hist[${bin}]=$((hist[${bin}] + 1))
    done
    
    # Normalizing factor
    sum=0
    for ((i=0; i<bins; i++)); do sum=$((sum + hist[$i])); done
    
    # Print vertical histogram
    # To keep the height dynamic we make so that the highest probability occupies 100 characters
    max_prob=$(printf '%s\n' "${hist[@]}" | sort -n | tail -1)
    max_prob=$(echo "${max_prob} / ${sum}" | bc -l )
    max_height=100

    # where we are going to print average and std dev
    ave_pos=$((max_height+30))
    ave_time=$(grep -m 1 "Average round trip time" ${output} | awk '{print $6}')
    std_time=$(grep -m 1 "Std dev round trip time" ${output} | awk '{print $7}')
    # Compute where average falls
    ave_bin=$(echo "scale=0; (${ave_time} - ${min}) / ${bin_width}" | bc)
    if [ ${ave_bin} -ge ${bins} ]; then ave_bin=$((bins - 1)); fi

    # Compute where stds fall
    std_bin_m=$(echo "scale=0; (${ave_time} - ${std_time} - ${min}) / ${bin_width}" | bc)
    if [ ${std_bin_m} -lt 0 ]; then std_bin_m=0; fi
    std_bin_p=$(echo "scale=0; (${ave_time} + ${std_time} - ${min}) / ${bin_width}" | bc)
    if [ ${std_bin_p} -ge ${bins} ]; then std_bin_p=$((bins - 1)); fi
    # print std deviations only if different from average
    if [[ ${std_bin_p} -ne ${ave_bin} && ${std_bin_m} -ne ${ave_bin} ]]; then
        do_std=1
    else
        do_std=0
    fi

    # Histogram header
    echo ""
    echo "|  Time  | Count | %Prob | Histogram"
    for ((i=0; i<bins; i++))
    do
        bin_count=${hist[${i}]}
        # Scale probability to max_height
        height=$(echo "${bin_count} * ${max_height} / ${max_prob} / ${sum}" | bc -l)
	# Convert to integer to always have integer number of characters 
        height=${height%.*}
        # Prepare string with bin centers
	bin_center=$(echo "${min}+${bin_width}*(${i} + 0.5)" | bc -l)
	str1=$(printf "| %6.3f |" ${bin_center})
	# Prepare string with counts
	str2=$(printf "%5g  |" ${bin_count})
	# Prepare string with probabilities
	str3=$(printf "  %4.1f |" $(echo "${bin_count} / ${sum} * 100" | bc -l))
	# Prepare string with histogram
	str4=$(printf " %s\n" "$(printf 'x%.0s' $(seq 1 ${height}))")
        # Full string
        full_string="${str1}${str2}${str3}${str4}"
	# report average
	if [ ${i} -eq ${ave_bin} ]; then
	    padded_string=$(printf "%-${ave_pos}s" "${full_string}")
	    full_string="${padded_string} < Average"
	fi
	# report std dev
	if { [[ ${i} -eq ${std_bin_m} ]] || [[ ${i} -eq ${std_bin_p} ]]; } && [[ ${do_std} -eq 1 ]]; then
            padded_string=$(printf "%-${ave_pos}s" "${full_string}")
            full_string="${padded_string} *"
        fi
        # report std dev bars
        if { [[ ${i} -gt ${std_bin_m} ]] && [[ ${i} -lt ${std_bin_p} ]] && [[ ${i} -ne ${ave_bin} ]] && [[ ${do_std} -eq 1 ]] }; then
            padded_string=$(printf "%-${ave_pos}s" "${full_string}")
            full_string="${padded_string} |"
        fi
	echo "${full_string}"
    done
fi

exit;
