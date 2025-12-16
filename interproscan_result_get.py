def find_proteins_with_multiple_domains(input_file, output_file):
    domain_dict = {}
    
    with open(input_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 5:
                continue
            protein_id = parts[0]
            domain_id = parts[4]
            
            if protein_id not in domain_dict:
                domain_dict[protein_id] = set()
            domain_dict[protein_id].add(domain_id)
    
    with open(output_file, 'w') as out_f:
        for protein_id, domains in domain_dict.items():
            if len(domains) > 1:
                out_f.write(protein_id + '\n')

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 3:
        print("Usage: python script.py input_file output_file")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    find_proteins_with_multiple_domains(input_file, output_file)
