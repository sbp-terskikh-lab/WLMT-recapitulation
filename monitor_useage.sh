#!/bin/bash

if [ $# -eq 0 ]
  then
    echo "No script specified. Usage: ./monitor_usage.sh <script>"
    exit 1
fi

# Start monitoring system resources
sar -u -r -d 1 > system_monitoring.log &

# Run the script
bash "$1"

# Stop monitoring system resources
killall sar

# Calculate max RAM, CPU, and disk space used
max_ram=$(grep -oP '\d+\.\d+' system_monitoring.log | awk '{print $1}' | sort -n | tail -1)
max_cpu=$(grep -oP '\d+\.\d+' system_monitoring.log | awk '{print $2}' | sort -n | tail -1)
# Calculate total disk space used
total_disk=$(grep -oP '\d+\.\d+' system_monitoring.log | awk '{print $3}' | paste -sd+ | bc)
total_disk=$(grep -oP '\d+\.\d+' system_monitoring.log | awk '{print $3}' | paste -sd+ | sed 's/^/0+/g' | bc)
# Save to log file
echo "Max RAM used: $max_ram" >> system_monitoring.log
echo "Max CPU used: $max_cpu" >> system_monitoring.log
echo "disk space used: $total_disk" >> system_monitoring.log

# Print
echo "Max RAM used: $max_ram"
echo "Max CPU used: $max_cpu"
echo "disk space used: $total_disk"

