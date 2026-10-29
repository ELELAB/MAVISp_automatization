import sys

marinara_list = sys.argv[1]

three_letter ={'V':'Val', 'I':'Ile', 'L':'Leu', 'E':'Glu', 'Q':'Gln', \
'D':'Asp', 'N':'Asn', 'H':'His', 'W':'Trp', 'F':'Phe', 'Y':'Tyr',    \
'R':'Arg', 'K':'Lys', 'S':'Ser', 'T':'Thr', 'M':'Met', 'A':'Ala',    \
'G':'Gly', 'P':'Pro', 'C':'Cys'}

def convert_code(x):
    out = ''
    for val in x:
        if val in three_letter.keys():
            val = val.replace(val, three_letter[val])
        else:
            val = val
        out += val
    return out

with open('input.csv', 'w') as m:
    m.write("name;site;type;function;reference" + '\n')
    with open(marinara_list, 'r') as vus:
        for row in vus:
            m.write(row.strip('\n') + ';p.'+ convert_code(row.rstrip('\n')) + ';mutation;;' + '\n')
