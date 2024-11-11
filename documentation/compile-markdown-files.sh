#!/bin/bash
#
# Use this script to compile all markedown files in this folder

for fname in "input-description.md" "configuration-options.md" "output-description.md"; do
    if [ $fname == "configuration-options.md" ]; then
        fontsize=10
    else
        fontsize=12
    fi
    echo "Compiling $fname with font size $fontsize ..."
    pandoc \
        -V 'fontfamily:dejavu' \
        -V 'mainfont:DejaVuSans' \
        -V 'sansfont:DejaVuSans' \
        -V 'monofont:DejaVuSansMono' \
        -V 'mathfont:TeXGyreDejaVuMath' \
        -V "fontsize=${fontsize}pt"\
        -V 'geometry:margin=25mm' \
        -V 'geometry:a4paper' \
        -V 'geometry:landscape' \
        -o ${fname::-2}pdf $fname
done
