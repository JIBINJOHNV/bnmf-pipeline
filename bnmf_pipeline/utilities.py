import os
import sys
import pandas as pd


def validate_trait_file(trait_file):
    """
    Validates trait file:
      - Can be read
      - Contains required columns
      - All trait_input_path files exist

    Returns:
        trait_file_df (pd.DataFrame)
    """

    # ----------------------------
    # Load file
    # ----------------------------
    try:
        trait_file_df = pd.read_csv(trait_file)
    except Exception as e:
        print(f"[ERROR] Failed to read trait file: {trait_file}")
        print(f"Reason: {e}")
        sys.exit(1)

    # ----------------------------
    # Check required columns
    # ----------------------------
    required_cols = {"trait_sample_id", "trait_input_path"}
    missing_cols = required_cols - set(trait_file_df.columns)

    if missing_cols:
        print(f"[ERROR] Missing required columns: {missing_cols}")
        sys.exit(1)

    # ----------------------------
    # Check input files exist
    # ----------------------------
    missing_files = [
        path for path in trait_file_df["trait_input_path"]
        if not os.path.exists(path)
    ]

    if missing_files:
        print("[WARNING] The following trait_input_path files do NOT exist:")
        for f in missing_files:
            print(f"   - {f}")
        sys.exit(1)

    print("[INFO] Trait file validation successful.")
    return trait_file_df