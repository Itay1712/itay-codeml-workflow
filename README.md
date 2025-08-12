This guide explains how to set up the environment and run the Codeml pipeline in test mode using a provided configuration script and Conda environment.

---

## Step 1: Set Up the Conda Environment

- Install Anaconda/Miniconda and create the `codeml_env` conda environment using `codeml_env.yaml`.
```bash
conda env create -f codeml_env.yaml \
&& conda activate codeml_env
```

---

## Step 2: Run the Pipeline in Test Mode

Once the environment is active, run the pipeline with the `config_test.yaml` test configuration:

```bash
f=codeml_pipeline.sh && sed -i 's/\r$//' "$f" && chmod +x "$f" && ./"$f" "config_test.yaml"
```

The configuration file now accepts an `OUTPUT_DIR` field defining the base folder where all pipeline results will be written. The pipeline automatically creates `processed` and `codeml` subdirectories inside this location.

---

## Expected Directory Structure

Your working directory should look like this:

```
project/
├── input/                     # Input that the pipeline uses
├── OUTPUT_DIR/                # Base directory for pipeline output (set in config)
│   ├── processed/             # Full output files from the pipeline
│   └── codeml/
│       ├── input/             # Codeml-specific input files selected from processed/
│       └── output/            # Codeml analysis results
├── codeml_pipeline.sh         # Main pipeline script
├── config_test.yaml           # Configuration file for test mode
└── codeml_env.yaml            # Conda environment definition
```

---

## Output Organization

Inside the `OUTPUT_DIR/codeml/output/` directory, results will be stored in a folder named after the `GROUP` variable from your config file.

For example, in a test run (that uses config_test.yaml) you will have:

```bash
GROUP="Astroviridae_43"
CODEML_RESULTS_DIR="${OUTPUT_DIR}/codeml/output/${GROUP}_test_$(date +%Y%m%d_%H%M%S)"
```

This creates a subfolder with the current time (to prevent overwrite) like:

```
OUTPUT_DIR/codeml/output/Astroviridae_43_test_20250614_111544/
├── Positive/      # Codeml branch-site analysis under positive selection
│   └── *.mlc      # Output files with likelihood model results
├── Null/          # Codeml null model for comparison
│   └── *.mlc
```

The `.mlc` files in these two folders are used to assess selection pressure using likelihood ratio tests.
