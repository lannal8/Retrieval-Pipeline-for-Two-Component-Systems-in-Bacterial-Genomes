import argparse

def main():
    parser = argparse.ArgumentParser(description='Process protein IDs based on domain start positions and c-Evalue.')
    parser.add_argument('-i', required=True, help='Input string for file names')
    args = parser.parse_args()
    a = args.i

    hatpase_list_file = f'08_HK_and_RR_list/{a}_HK.list'
    rr_list_file = f'08_HK_and_RR_list/{a}_RR.list'
    hatpase_out_file = f'03_hmmscan/{a}_HATPase.out'
    rr_out_file = f'03_hmmscan/{a}_RR.out'
    hhk_output_file = f'08_HK_and_RR_list/{a}_HHK.list'
    hrr_output_file = f'08_HK_and_RR_list/{a}_HRR.list'

    with open(hatpase_list_file, 'r') as f:
        hatpase_ids = set(line.strip() for line in f)
    with open(rr_list_file, 'r') as f:
        rr_ids = set(line.strip() for line in f)

    common_ids = hatpase_ids & rr_ids

    hatpase_data = {}
    with open(hatpase_out_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 18:
                continue
            protein_id = parts[3]  
            if protein_id in common_ids:
                try:
                    c_evalue = float(parts[11])  
                    start_pos = int(parts[17])  
                    
                    if protein_id not in hatpase_data or c_evalue < hatpase_data[protein_id]['c_evalue']:
                        hatpase_data[protein_id] = {
                            'c_evalue': c_evalue,
                            'start_pos': start_pos
                        }
                except (ValueError, IndexError):
                    continue

    rr_data = {}
    with open(rr_out_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 18:
                continue
            protein_id = parts[3]
            if protein_id in common_ids:
                try:
                    c_evalue = float(parts[11])
                    start_pos = int(parts[17])
                    
                    if protein_id not in rr_data or c_evalue < rr_data[protein_id]['c_evalue']:
                        rr_data[protein_id] = {
                            'c_evalue': c_evalue,
                            'start_pos': start_pos
                        }
                except (ValueError, IndexError):
                    continue

    with open(hhk_output_file, 'w') as hhk_f, open(hrr_output_file, 'w') as hrr_f:
        for pid in common_ids:
            if pid in hatpase_data and pid in rr_data:
                if hatpase_data[pid]['start_pos'] < rr_data[pid]['start_pos']:
                    hhk_f.write(pid + '\n')
                else:
                    hrr_f.write(pid + '\n')

if __name__ == '__main__':
    main()
