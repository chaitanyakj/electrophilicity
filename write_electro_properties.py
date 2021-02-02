"""
Write global and local electrophilicity properties.

Usage:

python write_electro_properties.py -f path/to/<filename>.log

path/to/filename.log should be supplied relative to the position of
write_electro_properties.py script.

Output (in the same folder as <filename>.log):

path/to/<filname>_local_electro.csv
path/to/<filename>_global_electro.txt

In the example I've used for reference, "<filename>" is "thiazole".
"""

import argparse
import pandas as pd


parser = argparse.ArgumentParser(description='Write global and local electrophilicity properties.')
parser.add_argument('-p','--path', help='Path to .log file', required=True)


def get_homo_lumo(path):
    """Get homo and lumo from a .log file."""
    with open(path, "r") as f:
        raw_text = f.read().split("\n")
        pop_indices = [
            i
            for i, x in enumerate(raw_text)
            if "Population analysis using the SCF density." in x
        ]
        raw_text = raw_text[pop_indices[-1] :]
        homo = float(
            [line for line in raw_text if "Alpha  occ. eigenvalues" in line][
                -1
            ].split()[-1]
        )
        "Alpha virt. eigenvalues"
        lumo = float(
            [line for line in raw_text if "Alpha virt. eigenvalues" in line][0].split()[
                4
            ]
        )
        return homo, lumo


def get_properties(homo, lumo):
    """Get properties from homo and lumo."""
    kcal_mol_homo = 27.212 * homo
    kcal_mol_lumo = 27.212 * lumo

    electronegativity = (-1 * (kcal_mol_homo + kcal_mol_lumo)) / 2
    hardness = kcal_mol_lumo - kcal_mol_homo
    softness = 1 / (2 * hardness)
    global_electrophilicity = electronegativity ** 2 / (2 * hardness)
    global_nucleophilicity = 1 / global_electrophilicity

    return (
        electronegativity,
        hardness,
        softness,
        global_electrophilicity,
        global_nucleophilicity,
    )


def write_global_electro(path):
    """Write global electrophilicity properties."""
    (electronegativity,
        hardness,
        softness,
        global_electrophilicity,
        global_nucleophilicity) = get_properties(*get_homo_lumo(path))
    written_text = (f"electronegativity; {electronegativity}\n"
                    f"hardness; {hardness}\n"
                    f"softness; {softness}\n"
                    f"global_electrophilicity; {global_electrophilicity}\n"
                    f"global_nucleophilicity; {global_nucleophilicity}")
    with open(f"{path[:-4]}_global_electro.txt", "w") as f:
        f.write(written_text)


def get_charge_table(in_path):
    """Get charge table."""
    with open(in_path, "r") as f:
        raw_text = f.read().split("\n")
        table_start_indices = [
            i for i, x in enumerate(raw_text) if "Mulliken atomic charges:" in x
        ]
        raw_text = raw_text[table_start_indices[-1] + 2 :]
        end_index = [i for i, x in enumerate(raw_text) if x.startswith(" Sum")][0]
        table_text = [line.split()[1:] for line in raw_text[:end_index]]
        table_text = [line for line in table_text if line[0] != "H"]
        element_dict = dict()
        for i in range(len(table_text)):
            table_text[i][1] = float(table_text[i][1])
            try:
                element_dict[table_text[i][0]] += 1
            except:
                element_dict[table_text[i][0]] = 1
            table_text[i][0] = f"{table_text[i][0]}{element_dict[table_text[i][0]]}"
        return table_text


def combine_tables(neutral_charge_table, negative_charge_table, positive_charge_table):
    """Combine properties from all charge tables into one."""
    overall_table = dict()
    for i in range(len(neutral_charge_table)):
        if (
            neutral_charge_table[i][0] == negative_charge_table[i][0]
            and neutral_charge_table[i][0] == positive_charge_table[i][0]
        ):
            overall_table[neutral_charge_table[i][0]] = [
                neutral_charge_table[i][1],
                positive_charge_table[i][1],
                negative_charge_table[i][1],
                neutral_charge_table[i][1] - positive_charge_table[i][1],
                negative_charge_table[i][1] - neutral_charge_table[i][1],
            ]
        else:
            raise Exception("Tables do not correspond.")
    return overall_table


def write_local_electro(path):
    """Write local electrophilicity properties."""
    neutral_charge_table = get_charge_table(path)
    negative_charge_table = get_charge_table(f"{path[:-4]}-1.log")
    positive_charge_table = get_charge_table(f"{path[:-4]}+1.log")
    overall_table = combine_tables(
        neutral_charge_table, negative_charge_table, positive_charge_table
    )

    (
        _,
        _,
        _,
        global_electrophilicity,
        global_nucleophilicity,
    ) = get_properties(*get_homo_lumo(path))
    for key in overall_table.keys():
        overall_table[key] += [
            overall_table[key][2] * global_nucleophilicity,
            overall_table[key][2] * global_electrophilicity,
        ]
    overall_table = pd.DataFrame.from_dict(overall_table).T.reset_index()
    overall_table.columns = [
        "Element",
        "Neutral",
        "(N-1)e'S=cation(+1)",
        "N+1=anion(-1)",
        "f-",
        "f+",
        "local_nucleophilicity",
        "local_elelectrorphilicity",
    ]
    overall_table.to_csv(f"{path[:-4]}_local_electro.txt", index=False)


if __name__ == '__main__':
    path = vars(parser.parse_args())["path"]
    write_global_electro(path)
    write_local_electro(path)