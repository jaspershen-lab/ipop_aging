# iPOP Aging Project

Welcome to the iPOP Aging repository! This repository contains code and documentation related to the iPOP Aging project. 

## About iPOP Aging Project

iPOP Aging is a research project focused on studying the effects of aging on the human body. Our goal is to understand the molecular nonlinear changes that occur in the human body as it ages, and help us better understand the aging process.

![](Figure_1.jpg)

## Repository Structure

- `1-code/`: This directory contains the code for the project.
- `2-data/`: This directory contains the data for the project.
- `3-data_analysis/`: This directory contains the data analysis for the project.
- `4-manuscript/`: This directory contains the manuscript for the project.
- `5-summary/`: This directory contains the summary for the project.

## Contact

If you have any questions or feedback, please feel free to reach out to us at [xiaotao.shen@outlook.com](xiaotao.shen@outlook.com) and [mpsnyder@stanford.edu](mpsnyder@stanford.edu).

## For authors only

How to ignore the files > 5M

### Step 1: Create a Script to Ignore Large Files

Create a script file (e.g., ignore_large_files.sh):

```bash
#!/bin/bash

# Find all files larger than 5MB and add them to .gitignore
find . -type f -size +5M | sed 's|^\./||' >> .gitignore

# Remove duplicate entries from .gitignore
sort -u -o .gitignore .gitignore

echo "All files larger than 5MB have been added to .gitignore."

```

### Step 2: Make the Script Executable

Make sure the script is executable:

```bash
chmod +x ignore_large_files.sh
```

### Step 3: Run the Script

Run the script to update your .`gitignore`:

```bash
./ignore_large_files.sh

```