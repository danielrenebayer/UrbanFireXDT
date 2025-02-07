#!/bin/bash

# This script compares the output of a program (stored in directories S0001 to S0011) against a verified reference output.
# It identifies files that differ, formats them for better readability, and presents them in vimdiff for manual inspection.
# The script pauses before opening each file, allowing the user to review the filename before proceeding.
# After comparison, temporary formatted files are deleted to keep the workspace clean.

# Check if the directory argument is provided
if [ -z "$1" ]; then
    echo "Usage: $0 <directory between S0001 and S0011>"
    exit 1
fi

dir="$1"

test_dir="test-output/$dir"
verified_dir="test-output-verified/$dir"

mapfile -t diff_files < <(diff -qr --exclude=build_and_run_info.txt --exclude=runtime-information.csv --exclude=ST1-AllCUs-ts.csv --exclude=ev-details.csv "$test_dir" "$verified_dir" | grep "differ" | sed -E 's/^Files (.*) and .* differ$/\1/')

# Loop through all differing files
for file in "${diff_files[@]}"; do
    relative_path="${file#$test_dir/}"
    verified_file="$verified_dir/$relative_path"

    tabbed_test_file="$file-tabbed"
    tabbed_verified_file="$verified_file-tabbed"

    # Print filename and wait for user input
    echo "Showing differences for: $relative_path"
    read -n 1 -s -r -p "Please press any key to continue..."
    echo ""

    # Convert CSV to tabbed format
    cat "$file" | column -t -s ',' > "$tabbed_test_file"
    cat "$verified_file" | column -t -s ',' > "$tabbed_verified_file"

    # Open in vimdiff
    vimdiff -o "$tabbed_test_file" "$tabbed_verified_file"

    # Clean up
    rm -f "$tabbed_test_file" "$tabbed_verified_file"
done

