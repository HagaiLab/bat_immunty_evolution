import os
import re
from scipy.stats import chi2
from statsmodels.stats.multitest import multipletests

def read_codeml_file(filename):
    """Reads in the contents of a codeML output file and extracts the Lrt value."""
    with open(filename, 'r') as f:
        contents = f.read()
        # Extract the lrt value from the output
        m = re.search(r"lnL\(.+?\):\s+([-0-9\.]+)", contents)
        if m:
            return float(m.group(1))
        else:
            print(f"No 'lnL' pattern match in file: {filename}")
    return None

def calculate_lrt(no_pos_file, pos_file):
    """Calculates the LRT between the two models."""
    lnL_no_pos = read_codeml_file(no_pos_file)
    lnL_pos = read_codeml_file(pos_file)
    if lnL_no_pos is None or lnL_pos is None:
        return None
    lrt = 2 * (lnL_pos - lnL_no_pos)
    return lrt

def main():
    directory = '/tzachi_storage/yanl/bat_genes/MAFFT_phyml/Tree/codeml/'
    
    # Files that pass the LRT test
    passed_file_path = '/tzachi_storage/yanl/bat_genes/MAFFT_positive_selection.txt'
    processed_prefixes = set()
    # Lists to store the p-values and LRT values
    p_values = []
    lrt_values = []
    
    for filename in os.listdir(directory):
        index = filename.find("_MSA")
        prefix = filename[:index]
        if filename.endswith('_MSA_no_positive.codeml'):
            # prefix has already been processed?
            if prefix not in processed_prefixes:
                processed_prefixes.add(prefix)
                # Check if there is a corresponding positive file with the same prefix
                pos_filename = prefix + '_MSA_positive.codeml'
                if os.path.exists(os.path.join(directory, pos_filename)):
                    lrt = calculate_lrt(os.path.join(directory, filename), os.path.join(directory, pos_filename))
                    if lrt is not None:
                        p_value = chi2.sf(lrt, df=1)
                        # Save the p-value and LRT value for FDR correction
                        p_values.append(p_value)
                        lrt_values.append((prefix, lrt))  # Store the prefix along with the LRT value
    
    # Apply FDR correction to the p-values
    reject, p_values_corrected, _, _ = multipletests(p_values, alpha=0.05, method='fdr_bh')
    
    with open(passed_file_path, 'w') as passed_file:
        # Iterate over the LRT values, corrected p-values, and prefixes
        for (prefix, lrt), p_value_corrected in zip(lrt_values, p_values_corrected):
            # Calculate the critical chi-square value for p=0.95 and df=1
            critical_value = chi2.ppf(0.95, df=1)
            # Check if LRT is higher than the critical chi-square value and the corrected p-value is significant
            if lrt > critical_value and p_value_corrected < 0.05:
                passed_file.write(f"{prefix}\t{lrt:.3f}\n")

if __name__ == '__main__':
    main()
