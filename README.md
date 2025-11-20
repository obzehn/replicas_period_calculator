# README
Takes log files from GROMACS for runs with replica exchange and computes the average times for replica zero to complete a round trip. The script produces some general statistics about the log from GROMACS, and an additional file that contains all the individual single round trips of replica zero for further analysis. This can be used in couple with `demux.pl` in the GROMACS scripts directory for more in depth analysis of replica exchange behaviour.

To use, just download the bash script and make it executable with `chmod +x`. It can be used in two ways:
1. Pass the log file to the script: `./replica_round_trip_extracter.sh <gromacs_log_file.log>`
2. Request an additional histogram with X bins on screen (mostly useless but cool): `./replica_round_trip_extracter.sh <gromacs_log_file.log> -histo X`
