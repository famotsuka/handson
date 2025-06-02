#Script to generate two plates with mixied mutagenesis primers (Fw and Rv)

import pandas as pd
import numpy as np
from Bio.Data import CodonTable


def OVP_FWprimers(FW_primers):
    """
    Function to generate a df with the Fw OVP primers
    """   
    # Filter the primers and empty columns
    FW_primers_filtered = FW_primers[FW_primers.iloc[:, 6].astype(str).str.contains("Fw_Frag_gg_", na=False)]
    # Further filter to keep only those matching the correct pattern: Fw_Frag_OVP_A_ or Fw_Frag_OVP_B_
    columns_to_drop = [
        "Payment Method", "Plate Barcode", "Sales Order #", "Reference #",
        "Measured Concentration uM", "Final Volume L",
        "Extinction Coefficient L/(mole·cm)", "Tm", "Well Barcode"
    ]
    FW_primers_filtered = FW_primers_filtered.drop(columns=[col for col in columns_to_drop if col in FW_primers_filtered.columns])

    return FW_primers_filtered

def OVP_RVprimers(RV_primers):
    """
    Function to generate a df with the Rv OVP primers
    """    
    # Filter the primers and empty columns
    RV_primers_filtered = RV_primers[RV_primers.iloc[:, 6].astype(str).str.contains("Rv_Frag_gg_", na=False)]
    # Further filter to keep only those matching the correct pattern: Fw_Frag_OVP_A_ or Fw_Frag_OVP_B_
    columns_to_drop = [
        "Payment Method", "Plate Barcode", "Sales Order #", "Reference #",
        "Measured Concentration uM", "Final Volume L",
        "Extinction Coefficient L/(mole·cm)", "Tm", "Well Barcode"
    ]
    RV_primers_filtered = RV_primers_filtered.drop(columns=[col for col in columns_to_drop if col in RV_primers_filtered.columns])

    return RV_primers_filtered

def generate_columnwise_wells():
    rows = ['A','B','C','D','E','F','G','H']
    cols = list(range(1, 13))  # 1 to 12
    return [f"{row}{col}" for col in cols for row in rows]


def create_transfer_dataframe(plate1,mutants):
    """
    Function to create a transfer dataframe for mixing primers
    """

    # Build reverse mapping and residue:codon dict. aa => [codons]
    standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
    aa_to_codons = {}
    for codon, aa in standard_table.forward_table.items():
        aa_to_codons.setdefault(aa, []).append(codon.upper())

    resulting_codons = {}
    for res_codons in mutants:
        aa = res_codons[-1]  # last character is the amino acid (e.g., 'R' in '5R')
        codons = aa_to_codons.get(aa[-1].upper(), [])
        resulting_codons[res_codons] = codons

    # Create the list of codons per residue:
    codon_list = []
    for mut in mutants:
        if mut in resulting_codons:
            codons = resulting_codons[mut]
            combined = [f"{mut}_{codon}" for codon in codons]
            codon_list.extend(combined)

    print("Codon combinations:")
    print(codon_list)   
    
    # Merge all the plates into one dataframe:
    FW_dfs_with_labels = [
    ("IDT_OVP", plate1),
    ("Fw_5inter", plate3),
    ("all15to41", plate5)
    ]
    RV_dfs_with_labels = [
    ("IDT_OVP", plate2),
    ("Rv_5inter", plate4),
    ("all15to41", plate6)
    ]

    all_fw_df = pd.concat(
        [df.assign(**{"Source Plate": label}) for label, df in FW_dfs_with_labels],
        ignore_index=True
    )
    all_rv_df = pd.concat(
        [df.assign(**{"Source Plate": label}) for label, df in RV_dfs_with_labels],
        ignore_index=True
    )
 
    # Create the transfer dataframe
    vol = 5
    Destination_plate = "SPreaction"
    columns = ['Name', 'Source Plate', 'Source Well', 'Destination Plate', 'Destination Well', 'Volume']
    fw_data = []
    rv_data = []
    destination_wells = generate_columnwise_wells()
    well_index = 0
    well_index_Rv = 0
    for codon_id in codon_list:
        matches = all_fw_df[all_fw_df["Name"].str.contains(codon_id)]
        if not matches.empty:
            for _, row in matches.iterrows():
                fw_data.append({
                    "Name": row["Name"],
                    "Source Plate": row["Source Plate"],
                    "Source Well": row["Well Position"],
                    "Destination Plate": Destination_plate,
                    "Destination Well": destination_wells[well_index],
                    "Volume": vol
                })
                well_index += 1

    for codon_id in codon_list:
        matches = all_rv_df[all_rv_df["Name"].str.contains(codon_id)]
        if not matches.empty:
            for _, row in matches.iterrows():
                rv_data.append({
                    "Name": row["Name"],
                    "Source Plate": row["Source Plate"],
                    "Source Well": row["Well Position"],
                    "Destination Plate": Destination_plate,
                    "Destination Well": destination_wells[well_index_Rv],
                    "Volume": row["nmoles"] * 10
                })
                well_index_Rv += 1
    
    print(type(fw_data))
    DF_final = pd.concat([
    pd.DataFrame(fw_data, columns=columns),
    pd.DataFrame(rv_data, columns=columns)
    ], ignore_index=True)
    
    return DF_final


def dissolve_water(plate_df):
    """
    Create a transfer plan for dissolving primers with water.
    
    Each primer will receive nmoles * 10 µL of water.
    The water is taken from a column-wise ordered 96-well plate.
    """
    destination_wells = generate_columnwise_wells()
    Destination_plate = "stock_primers"

    transfer_rows = []
    well_index = 0

    for _, row in plate_df.iterrows():
        transfer_rows.append({
            "Name": row["Sequence Name"],
            "Source Plate": "Water",
            "Source Well": destination_wells[well_index],
            "Destination Plate": Destination_plate,
            "Destination Well": row["Well Position"],
            "Volume": round(row["nmoles"] * 10, 1)
        })
        well_index += 1

    dissolve_df = pd.DataFrame(transfer_rows)

    
    return dissolve_df

def create_working_concentration(plate_df, conc, vol_final):
    init_conc = 100  # initial concentration in µM
    destination_plate = "working_primers"
    destination_wells = generate_columnwise_wells()
    
    working_rows = []
    well_index = 0

    for _, row in plate_df.iterrows():
        working_rows.append({
            "Name": row["Name"],
            "Source Plate": row["Destination Plate"],  # source is previous destination
            "Source Well": row["Destination Well"],
            "Destination Plate": destination_plate,
            "Destination Well": destination_wells[well_index],
            "Volume": round((conc * vol_final) / init_conc, 1)  # dilution calc
        })
        well_index += 1

    working_df = pd.DataFrame(working_rows)
    return working_df    


import pandas as pd

def combine_fw_rv_primers(df, destination_plate):
    """
    Combines Rv and Fw primers into a new DataFrame with specified destination wells and volume.
    
    Parameters:
        input_csv (str): Path to the CSV file with primer info.
        destination_plate (str): Name of the destination plate (e.g. 'comb_frags_gg').
        destination_wells (list): List of wells for the combined output (e.g. ['A1', 'B1', ...]).
    
    Returns:
        pd.DataFrame: A new DataFrame with combined primer information.
    """
    # Read the CSV
    # Split into forward and reverse primer rows
    fw_df = df[df['Name'].str.contains("Fw_Frag_gg_2_", na=False)]
    rv_df = df[df['Name'].str.contains("Rv_Frag_gg_1_", na=False)]

    transfers = []
    well_index = 0
    destination_wells = generate_columnwise_wells()

    for _, rv_row in rv_df.iterrows():
        for _, fw_row in fw_df.iterrows():
            if well_index >= len(destination_wells):
                raise ValueError("Not enough destination wells for all combinations.")

            destination_well = destination_wells[well_index]
            combined_name = f"{fw_row['Name']}__{rv_row['Name']}"

            # Add transfer from Fw primer
            transfers.append({
                "Name": combined_name,
                "Source Plate": fw_row["Destination Plate"],
                "Source Well": fw_row["Destination Well"],
                "Destination Plate": destination_plate,
                "Destination Well": destination_well,
                "Volume": 5
            })

            # Add transfer from Rv primer
            transfers.append({
                "Name": combined_name,
                "Source Plate": rv_row["Destination Plate"],
                "Source Well": rv_row["Destination Well"],
                "Destination Plate": destination_plate,
                "Destination Well": destination_well,
                "Volume": 5
            })

            well_index += 1

    return pd.DataFrame(transfers)

