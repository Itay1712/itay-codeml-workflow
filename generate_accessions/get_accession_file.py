import pandas as pd
import os, sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
DATASET_DIR = os.path.join(PROJECT_ROOT, "generate_accessions")

def extract_values(input_file, sheet_name, column1, value1, columns, header=0):
    df = pd.read_excel(input_file, sheet_name=sheet_name, header=header)
    filtered_df = df[df[column1] == value1]
    result_values = filtered_df[columns].to_dict(orient='records')
    return result_values

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(f"Usage: {sys.argv[0]} <group_name> <output_file>")
    
    group=sys.argv[1]
    output_file=sys.argv[2]
    input_file = os.path.join(DATASET_DIR, "41559_2024_2353_MOESM5_ESM.xlsx")

    ## step 1: extract accession + host genus from Supplementary Table 1
    values_table1 = extract_values(
        input_file=input_file,
        sheet_name="Supplementary Table 1",
        column1="Viral clique",
        value1=group,
        columns=["Accession", "Host genus"],
        header=1
    )
    df1 = pd.DataFrame(values_table1)
    df1 = df1.rename(columns={"Host genus": "Host"})

    ## step 2: extract accession + ancestral node from Supplementary Table 2
    values_table2 = extract_values(
        input_file=input_file,
        sheet_name="Supplementary Table 2",
        column1="Clique name",
        value1=group,
        columns=["Tip name", "Ancestral node"],
        header=4
    )
    df2 = pd.DataFrame(values_table2)
    df2 = df2.rename(columns={"Tip name": "Accession"})
    
    ## step 3: merge by accession
    merged_df = pd.merge(df1, df2, on="Accession", how="left")

    ## step 4: fill missing ancestral node with "Root"
    merged_df["Ancestral node"] = merged_df["Ancestral node"].fillna("Root")

    ## step 5: remove rows where host is missing or empty
    for _, row in merged_df.iterrows():
        host = str(row["Host"]).strip()
        if not host or host.lower() == "nan":
            print(f"!!!  Dropped due to missing Host: {row['Accession']}")

    merged_df = merged_df[merged_df["Host"].notna()]
    merged_df = merged_df[merged_df["Host"].astype(str).str.strip() != ""]
    merged_df = merged_df.reset_index(drop=True)

    ## step 6: reorder columns and export
    final_df = merged_df[["Ancestral node", "Accession", "Host"]]
    final_df.to_csv(output_file, sep="\t", index=False)

    print("Accessions file output:\n", final_df)
