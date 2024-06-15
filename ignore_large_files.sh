#!/bin/bash

# Find all files larger than 5MB and add them to .gitignore
find . -type f -size +5M | sed 's|^\./||' >> .gitignore

# Remove duplicate entries from .gitignore
sort -u -o .gitignore .gitignore

echo "All files larger than 5MB have been added to .gitignore."
