import mix_primers
import pandas as pd
import re
import os

# list of residues to mutate with combined fragments
muts = ['40L','52S']

def clean_well_position(well):
    match = re.match(r"([A-H])0([1-9])$", well)
    if match:
        return f"{match.group(1)}{match.group(2)}"
    return well  # Return unchanged if already correct


def main():
    print("\n - Hello from double-syn-muts-ddla :D")
    idt_primers = pd.read_csv("IDT_Plate_Specs.csv")
    idt_primers["Well Position"] = idt_primers["Well Position"].apply(clean_well_position)
    print("\nloaded IDT primers in a plate:")
    print(idt_primers.head())

    os.makedirs("outputs", exist_ok=True)

    result_Fw_IDT_OVP = mix_primers.OVP_FWprimers(idt_primers)
    result_Rv_IDT_OVP = mix_primers.OVP_RVprimers(idt_primers)
    #print(result_Fw_IDT_OVP.head())
    #print(result_Rv_IDT_OVP.head())

    # Prepare Stock solution of primers:
    DF_stock = pd.concat([result_Fw_IDT_OVP, result_Rv_IDT_OVP], ignore_index=True)
    Stock_primers = mix_primers.dissolve_water(DF_stock)
    print("\n - Mix Primers. Input to make stock primers of 100uM in the original plate using Biomek i7:")
    print(Stock_primers.head())
    # Save stock primers to CSV files:
    output_file1="outputs/stock_primers.csv"
    output_file2="outputs/stock_primers_biomek.csv"
    Stock_primers.to_csv(output_file1, index=False)
    dissolve_df_biomek = Stock_primers.drop(columns=["Name"])
    dissolve_df_biomek.to_csv(output_file2, index=False)
    print(f"\nFiles available in: {output_file1} and {output_file2}")

    # Prepare working concentration of Fw and Rv primers in uM:
    working_primers = mix_primers.create_working_concentration(Stock_primers, conc=10, vol_final=200)
    print("\n - Working primers. Input to make working solution of primers (10 uM) in microtitter plate using Biomek i7:")
    print(working_primers.head())
    # Save stock primers to CSV files:
    output_file3 = "outputs/working_primers.csv"
    output_file4 = "outputs/working_primers_biomek.csv"
    working_primers.to_csv(output_file3, index=False)
    working_primers_biomek = working_primers.drop(columns=["Name"])
    working_primers_biomek.to_csv(output_file4, index=False)
    print(f"\nFiles available in: {output_file3} and {output_file4}")


    combine_frags = mix_primers.combine_fw_rv_primers(working_primers, destination_plate="comb_frags_gg")
    print(f"\n - Combine fragments. Input to combine 2 fragments for gene assemble using Biomeki7:")
    print(combine_frags.head())
    output_file5 = "outputs/combined_frags.csv"
    output_file6 = "outputs/combined_frags_biomek.csv"
    combine_frags.to_csv(output_file5, index=False)
    combine_frags_biomek = combine_frags.drop(columns=["Name"])
    combine_frags_biomek.to_csv(output_file6, index=False)
    print(f"\nFiles available in: {output_file5} and {output_file6}\n")


    print(r"Congratz \o/ you have successfully generated the primers and the input files necessary \n" \
    r"to make double mutants in your gene using the Biomek i7 :D")


if __name__ == "__main__":
    main()

