"""
Write global and local electrophilicity properties from NBO Charges.

Usage:

python write_electro_nbo.py -p path/to/<filename>.log

path/to/filename.log should be supplied relative to the position of
write_electro_properties.py script.

Output (in the same folder as <filename>.log):

path/to/<filname>_electro.csv

In the example I've used for reference, "<filename>" is "thiazole".
"""

import argparse
import pandas as pd


# Constants to help parse charge table.
SECTION_HEADER = "******************************Gaussian NBO Version 3.1******************************"
TABLE_HEADER = " Summary of Natural Population Analysis:"
TABLE_END = " ======================================================================="


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
	written_text = (f"{path[:-4]}\n"
                    f"electronegativity,{electronegativity:.4f}\n"
					f"hardness,{hardness:.4f}\n"
					f"softness,{softness:.4f}\n"
					f"global_electrophilicity,{global_electrophilicity:.4f}\n"
					f"global_nucleophilicity,{global_nucleophilicity:.4f}"
                    "\n\n")
	with open(f"{path[:-4]}_electro.csv", "w") as f:
		f.write(written_text)


def get_charge_table(in_path):
	"""Get charge table."""
	with open(in_path, "r") as f:
		raw_text = f.read().split("\n")
		section_start_indices = [
			i for i, x in enumerate(raw_text) if SECTION_HEADER in x
		]
		raw_text = raw_text[section_start_indices[-1] :]
		table_start_indices = [i for i, x in enumerate(raw_text) if TABLE_HEADER in x]

		raw_text = raw_text[table_start_indices[0] + 6 :]
		end_index = [i for i, x in enumerate(raw_text) if x.startswith(TABLE_END)][0]
		raw_text = raw_text[:end_index]

		table_text = [line.split() for line in raw_text[:end_index]]
		table_text = [
			[f"{entry[0]}{entry[1]}", float(entry[2])]
			for entry in table_text
			if entry[0] != "H"
		]
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
	overall_table = pd.DataFrame.from_dict(overall_table).T.reset_index().round(4)
	overall_table.columns = [
		"Element",
		"Neutral",
		"(N-1)e'S=cation(+1)",
		"N+1=anion(-1)",
		"f-",
		"f+",
		"local_nucleophilicity",
		"local_electrophilicity",
	]
	overall_table.to_csv(f"{path[:-4]}_electro.csv", index=False, mode="a")


if __name__ == '__main__':
	path = vars(parser.parse_args())["path"]
	write_global_electro(path)
	write_local_electro(path)