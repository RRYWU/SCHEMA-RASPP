import sys

def read_clustal(file_name):
    seqs = {}
    order = []
    with open(file_name, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('CLUSTAL'):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            name, seq = parts[0], parts[1]
            if name not in seqs:
                seqs[name] = ''
                order.append(name)
            seqs[name] += seq
    return seqs, order

def make_chimera(msa_dict, seq_order, chimera_code, xo_points):
    seq_length = len(msa_dict[seq_order[0]])
    xo_points = [0] + xo_points + [seq_length]

    chimera_seq = ''
    for i, parent_idx_char in enumerate(chimera_code):
        if not parent_idx_char.isdigit():
            continue  
        parent_idx = int(parent_idx_char) - 1  
        start = xo_points[i]
        end = xo_points[i + 1]
        fragment = msa_dict[seq_order[parent_idx]][start:end]
        chimera_seq += fragment
    return chimera_seq

def main():
    if len(sys.argv) < 5:
        print("Usage: python chimeraseq.py parent_msa.aln chimera_codes.txt xo_points.txt output.fasta")
        sys.exit(1)
    
    msa_file = sys.argv[1]
    chimera_file = sys.argv[2]
    xo_file = sys.argv[3]
    out_file = sys.argv[4]

    msa_dict, seq_order = read_clustal(msa_file)

    xo_points = []
    with open(xo_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            xo_points = [int(x) for x in line.split()]
            break

    chimera_codes = []
    with open(chimera_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                chimera_codes.append(line)

    with open(out_file, 'w') as f:
        for code in chimera_codes:
            chimera_seq = make_chimera(msa_dict, seq_order, code, xo_points)
            f.write(">chimera_%s\n%s\n" % (code, chimera_seq))

    print("{} chimera sequences written to {}".format(len(chimera_codes), out_file))
    print("Parent order used: {}".format(seq_order))

if __name__ == "__main__":
    main()
