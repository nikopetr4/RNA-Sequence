# Initializations
#rna_seq = "AUGGCCAUUGUAAUGGGCCGCUGA"
rna_seq = 'ACGUCGAUUCGAGCGAAUCGUAACGAUACGAGCAUAGCGGCUAGAC'
#rna_seq = 'ACAUGAUGGCCAUGU'
#rna_seq = 'ACCGGUAGU'
#rna_seq = 'AUGUGGCCAU'
#rna_seq = 'AAGACUUCGGAUCUGGCGACACCC'


def print_contents(opt_val):
    n = len(opt_val)
    for i in range(n-5-1, -1, -1):
        print(i, end=': ')
        for j in range(5, n):
            print(str(opt_val[i][j]).rjust(2), end=' ')
        print()
    print('j:', end=' ')
    for j in range(5, n):
        print(str(j).rjust(2), end=' ')
    print()
    print()


sequence = rna_seq
n = len(sequence)
opt_val = [[-1 for j in range(n)] for i in range(n)]
opt_val_pairs = [[[] for j in range(n)] for i in range(n)]

for i in range(n):
    for j in range(n):
        if i >= j - 4:
            opt_val[i][j] = 0

setAU = set('A')
setAU.add('U')
setCG = set('C')
setCG.add('G')
valid_pairs = [setAU, setCG]


def can_pair(base_1: str, base_2: str) -> bool:
    pair = set(base_1)
    pair.add(base_2)
    return pair in valid_pairs


def get_opt_val(sequence:str, n:int, opt_val, opt_val_pairs, i:int, j:int) -> int:
    if opt_val[i][j] == -1:
        if i >= j - 4:
            opt_val[i][j] = 0       # no pair can be made
        else:
            term2 = -1
            index = -1
            term1 = get_opt_val(sequence, n, opt_val, opt_val_pairs, i, j - 1)
            for t in range(i, j - 4):
                if can_pair(sequence[t], sequence[j]):
                    if t < 1:
                        val = 1 + get_opt_val(sequence, n, opt_val, opt_val_pairs, t + 1, j - 1)
                        if val > term2:
                            term2 = val
                            index = t
                    else:
                        val = 1 + get_opt_val(sequence, n, opt_val, opt_val_pairs, i, t - 1) + get_opt_val(sequence, n, opt_val, opt_val_pairs, t + 1, j - 1)
                        if val > term2:
                            term2 = val
                            index = t

            opt_val[i][j] = max(term1, term2)

            if term1 >= term2:
                opt_val_pairs[i][j].extend(opt_val_pairs[i][j - 1])
            else:
                opt_val_pairs[i][j].append([index, j])
                opt_val_pairs[i][j].extend(opt_val_pairs[index + 1][j - 1])
                opt_val_pairs[i][j].extend(opt_val_pairs[i][index - 1])
    return opt_val[i][j]


def get_dot_bracket_notation(n: int, sorted_pairs)->str:
    s = '.'*n
    slist = list(s)
    for pair in sorted_pairs:
        slist[pair[0]] = '('
        slist[pair[1]] = ')'
    print(slist)
    finString = "".join(slist)
    return finString

def compute_rna_secondary_structure(sequence):
    # Initializations
    n = len(sequence)
    opt_val = [[-1 for j in range(n)] for i in range(n)]
    opt_val_pairs = [[[] for j in range(n)] for i in range(n)]

    for i in range(n):
          for j in range(n):
            if i >= j - 4:
                opt_val[i][j] = 0


    for k in range(5, n):
          for i in range(0, n - k):
            j = i + k
            get_opt_val(sequence,n,opt_val,opt_val_pairs,i, j)
            print_contents(opt_val)

    for i in range(n):
        for j in range(n):
            print(opt_val_pairs[i][j])
    print(opt_val_pairs[0][n-1])

    sorted_pairs = sorted(opt_val_pairs[0][n-1])
    sorted_pair_bases = [(sequence[i], sequence[j]) for i, j in sorted_pairs]
    dot_bracket_notation = get_dot_bracket_notation(n , sorted_pairs)
    return dot_bracket_notation, sorted_pairs, sorted_pair_bases, opt_val, opt_val_pairs

compute_rna_secondary_structure(rna_seq)
