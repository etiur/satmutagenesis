def map_atom_string(atom_string, initial_pdb, prep_pdb):
    # read in user input
    with open(initial_pdb, "r") as initial:
        initial_lines = initial.readlines()

    # read in preprocessed input
    with open(prep_pdb, "r") as prep:
        prep_lines = prep.readlines()

    # split the atom string or retrieve from the number
    if len(atom_string.split(":")) == 3:
        chain, resnum, atom_name = atom_string.split(":")
        # extract coordinates from user input
        for i in initial_lines:
            if (i.startswith("HETATM") or i.startswith("ATOM")) and i[21].strip() == chain.strip() and i[
                    22:26].strip() == resnum.strip() and i[12:16].strip() == atom_name.strip():

                coords = i[30:54].split()

                # extract coordinates from preprocessed file
                for p in prep_lines:
                    if p[30:54].split() == coords:
                        new_atom_name = p[12:16].strip()
                        new_resnum = p[22:26].strip()
                        new_chain = p[21].strip()
                        after = "{}:{}:{}".format(new_chain, new_resnum, new_atom_name)
                        break
                break

    else:
        chain, resnum = atom_string.split(":")
        for i in initial_lines:
            if (i.startswith("HETATM") or i.startswith("ATOM")) and i[21].strip() == chain.strip() and i[
                                                                    22:26].strip() == resnum.strip():
                coords = i[30:54].split()

                # extract coordinates from preprocessed file
                for p in prep_lines:
                    if p[30:54].split() == coords:
                        new_resnum = p[22:26].strip()
                        new_chain = p[21].strip()

                        after = "{}:{}".format(new_chain, new_resnum)

                        break
                break

    return after



