# This folder helps us to preprocess the data of smiles

1) If needed to calculate the total number of heavy atoms in smiles : "atom_count.py" --->> python atom_count.py --input input_file --output output_file(in either.dat or .csv format)

2) If needed to calculate the logp values of smiles : "logp.py" --->> python atom_count.py --input input_file --output output_file(in either.dat or .csv format)

3) If you need to calculate all the values of smiles say for atom count less than particular number and what are the properties of that smiles : data_preprocess.py ------>> python data_preprocess.py --input input_file --output output_file(in either.dat or .csv format)

4) Once you done with analysis and if you have the desire output with say logP values and MR values then you can find out what the data stays in the range : "actual_mol.py" -----> python actual_mol.py --input input_file --data1 min_logp_value --data2 max_logp_value --data3 min_mr_value --data4 max_mr_value
