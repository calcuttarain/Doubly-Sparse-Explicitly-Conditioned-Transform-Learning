---

Run `./run_benchmark.m`. It outputs a folder with results containing a timestamp.

Run `./analyze_results.m` by specifying the parameters `timestamp`, `lambda_idx`, and `parameters_settings_idx` at the beginning of the script. It outputs a timestamped folder containing the generated plots, and displays an error report for any failed configurations.

---

Datasets:
1. **The original authors' dataset**, located in the `./data` folder.
2. **Partial DIV2K dataset**, located in the `./DIV2K_valid_HR_data` folder. This contains the first 10 images from the DIV2K High-Resolution validation set, grayscaled and center-cropped to 1024x1024. (The exact preprocessing pipeline is available in `./utils/download_preprocess_div2k_dataset.m`).

---

The pipeline tests the following configurations:
- both datasets
- various input values ($\lambda$) for the Bresler methods
- various representation (`T0`) and transform sparsity levels (`T1`)
- multiple patch sizes (and implicitly, transform sizes) (`n`)
- the Unstructured and Structured Conditioned methods against three Bresler methods: Unstructured Closed Form, Structured, and Structured Closed Form.

The exact values for these settings are defined at the beginning of `run_benchmark.m`.

---

**Note:** 
- This code uses `parfor` for parallel execution, which requires Parallel Computing Toolbox.
- Doubly Transform and Structured Transform terms are used interchangeably.

---