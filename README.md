# **FBA Model Conversion and Testing Repository**

This repository is a **test framework** designed to evaluate the conversion of metabolic models from MATLAB `.mat` format to a standardized `.txt` format for Flux Balance Analysis (FBA). The repository also includes functionality for setting custom flux bounds and applying parsimonious FBA (pFBA) to the models.

---

## **Features**

1. **MATLAB to TXT Conversion**:
   - Converts `.mat` models into `.txt` format, ensuring compatibility with downstream FBA workflows.
   - Provides a structured output for easier analysis and reproducibility.

2. **Flux Bound Customization**:
   - Allows setting specific upper bounds for **exchange reactions (EX_)** in the models.
   - Provides the flexibility to define bounds for both **FBA-specific reactions** and general metabolic reactions.

3. **pFBA Support**:
   - Offers the ability to enable parsimonious FBA (pFBA) for optimizing models:
     - **`geneOption = 0`**: Minimizes all reactions.
     - **`geneOption = 1`**: Minimizes only reactions associated with genes.
     - **`geneOption = 2`**: Minimizes only reactions not associated with genes.
     - **Default = -1**: pFBA is disabled for the model.

---

## **Repository Structure**

- **`/input/`**: Contains the `.mat` files for each metabolic model categorized by organism type.
- **`/epimod_FBAfunctions/`**: Includes R scripts for:
  - Model generation (`FBAgreatmodeClass.R`, `class_generation.R`).
  - Reading `.mat` files (`readMat.R`).
  - Setting bounds for exchange reactions (`ex_bounds_module.R`).
- **`/results/`**: Stores the converted `.txt` files and processed outputs.
- **`main.R`**: The primary script orchestrating the workflow.

---

## **Main Script (`main.R`)**

The `main.R` script handles the following tasks:

1. **Model Conversion**:
   - Loads models from `.mat` files using the `FBA4Greatmod.generation` function.
   - Configures biomass parameters (maximum, mean, and minimum).
   - Converts the models into `.txt` format for further analysis.

2. **File Management**:
   - Automatically moves `.rds` files generated during conversion to the appropriate directory for storage.

3. **Setting Upper Bounds**:
   - Extracts **exchange reactions (EX_)** from the `.txt` models.
   - Configures upper bounds for:
     - Non-FBA reactions (`non_fba_base_bound / count`).
     - FBA-specific reactions (`fba_upper_bound`).

4. **Optional pFBA Application**:
   - Enables parsimonious FBA by setting the `geneOption` flag, if specified.

5. **Execution Tracking**:
   - Records and reports the time taken to process each model and the total execution time.

---

## **How to Use**

1. **Setup**:
   - Place your `.mat` files in the `/input/<type>/` directory (e.g., `/input/ecoli/`).
   - Ensure the working directory is correctly set in the script.

2. **Run the Main Script**:
   - Execute `main.R` to process all specified models.
   - Monitor the output in the `/results/` directory.

3. **Customize Bounds**:
   - Modify the `run_full_ex_bounds` parameters in `main.R` to adjust reaction bounds as needed.

4. **Enable pFBA**:
   - Uncomment and set the `setPFbaGeneOption` function in `main.R` to enable pFBA with the desired configuration.

---

## **Contact**

For questions or feedback, feel free to reach out to the repository maintainer.

