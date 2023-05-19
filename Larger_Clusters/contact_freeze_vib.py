
with open('CONTCAR', 'r') as infile:
    line = infile.readlines()
#print(line)
#for i in range(len(line)):
#    print(i, line[i])
total_atom_line = line[6].strip()
atom_type_line = line[5].strip()
atom_count = total_atom_line.split()
atom_types = atom_type_line.split()
atom_count_int = list(map(int, atom_count))
#print(atom_count_int)
#print(atom_types)
tot_atoms = sum(atom_count_int)
tot_metal_atoms = atom_count_int[0]
tot_adsorb_atoms = tot_atoms - atom_count_int[0]
#print("tot_atoms", tot_atoms)
#print("tot_metal_atoms", tot_metal_atoms)
#print("tot_adsorb_atoms", tot_adsorb_atoms)

with open('CONTCAR_Freq_xyz', 'w') as outfile:
    for i in range(7):
        outfile.write(line[i])
    outfile.write("Selective dynamics \n")
    outfile.write(line[7])



sum_lines_metals = (8, 8+tot_metal_atoms)
sum_lines_adsorbs = (sum_lines_metals[1], sum_lines_metals[1] + tot_adsorb_atoms)
#print( sum_lines_metals, sum_lines_adsorbs)

for i in range(sum_lines_metals[0], sum_lines_metals[1]):
        with open('CONTCAR_Freq_xyz', 'a+') as outfile:
            new_line = line[i].rstrip() + str('  F F F') + '\n'
            outfile.write(new_line)

for i in range(sum_lines_adsorbs[0], sum_lines_adsorbs[1]):
        with open('CONTCAR_Freq_xyz', 'a+') as outfile:
            new_line = line[i].rstrip() + str('  T T T') + '\n'
            outfile.write(new_line)



