Here's a refined version of your GitHub README file with improved structure and formatting:

```md
# EnzBuilder: Simultaneously Improving Metabolic Model Reconstruction and Enzyme Function Annotation with Machine Learning Methods

EnzBuilder integrates predicted protein annotations directly into the metabolic reconstruction and refinement processes.

---

## Table of Contents
1. [About](#about)
2. [Features](#features)
3. [Installation](#installation)
4. [Usage](#usage)
5. [CLEAN Prediction](#clean-prediction)
6. [Contributing](#contributing)
7. [License](#license)
8. [Acknowledgments](#acknowledgments)

---

## About
EnzBuilder is a tool designed to bridge the gap between protein annotation and metabolic model reconstruction by leveraging machine learning methods.  
- Incorporates predicted enzyme functions into draft models.
- Enables iterative model refinement using AI predictions.
- Enhances accuracy and coverage of metabolic networks.

---

## Features
- **Seamless integration** of machine learning-based protein annotations.
- **Automated model gap-filling** based on functional predictions.
- **Compatibility** with common COBRA and SBML workflows.
- **Iterative refinement** of metabolic models to improve network coverage.

---

## Installation
To install and set up EnzBuilder, follow these steps:

```bash
# Clone the repository
git clone https://github.com/yourusername/EnzBuilder.git

# Navigate to the project directory
cd EnzBuilder

# Install dependencies
conda create --name enzbuilder --file enzbuilder.yml

# Activate the conda environment
conda activate enzbuilder
```

### Data Download
1. Download the required data file from [Google Drive](https://drive.google.com/file/d/1GFBLrw4uxEg1Ht67DdZn7K0g1hw9p8QN/view?usp=drive_link).
2. Extract the data by running:

```bash
cd data
tar -xzvf data.tar.gz
```

---

## Usage
### Basic Workflow
Run the following command to process your input file:

```bash
python funcarve_main.py \
  --input_file ../data/test/CP000148.1_t4_maxsep_df.pkl \
  --gram negative \
  --block_flage 0 \
  --flux_flage 0 \
  --name iLJ478_allec \
  --media default \
  --reward 0.1 \
  --iter 3 \
  --cpu 8 \
  --threshold 5 \
  --upper 15 \
  --lower 5 \
  --maxweight 100 \
  --minweight 0.0
```

---

## CLEAN Prediction
EnzBuilder relies on **CLEAN** for enzyme function predictions. Follow these steps to set up CLEAN:

1. Install CLEAN by following the instructions at [CLEAN GitHub](https://github.com/tttianhao/CLEAN).
2. Replace `CLEAN/app/src/CLEAN/infer.py` with the version in `src/infer.py` from EnzBuilder.
3. When running CLEAN, an additional `.pkl` file will be generated. Use this file as input for EnzBuilder.

---

## Contributing
We welcome contributions! If you'd like to contribute, please fork the repository, create a branch, and submit a pull request.

For bug reports or feature requests, please open an issue on the GitHub repository.

---

<!-- ## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details. -->

---

## Acknowledgments
- Thank you to all contributors and researchers who have helped shape EnzBuilder.
- Special thanks to the authors of CLEAN for their enzyme function prediction framework.
