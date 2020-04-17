
out_fn = "/Users/kmayerbl/active/cagalog/taxdump/taxid_to_sciname.csv"
fn = "/Users/kmayerbl/active/cagalog/taxdump/names.dmp"
with open(out_fn ,"w") as oh:
    oh.write(f'taxid,scientific_name\n')
    with open(fn, "r") as fh:
        c = 0
        for line in fh:
            line = line.strip().split("|")
            line = [x.replace("\t","") for x in line]
            if line[3] != "scientific name":
                continue
            oh.write(f'{line[0]},{line[1]}\n')

out_fn = "/Users/kmayerbl/active/cagalog/taxdump/taxid_to_sciname.csv"
taxid_to_sciname = dict()
with open(out_fn, "r") as fh:
    for line in fh:
        k,v = line.strip().split(",")
        taxid_to_sciname[k] = v







