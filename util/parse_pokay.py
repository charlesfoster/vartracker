import os
import pandas as pd
from urlextract import URLExtract
import glob
import sys

usage = "Usage: python parse_pokay.py pokay_data_dir output_file"

# Mapping of product names to nsp names
PRODUCT_TO_NSP = {
    "PLpro": "nsp3",
    "3CLpro": "nsp5",
    "3CL": "nsp5",
    "RdRp": "nsp12",
    "Hel": "nsp13",
    "ExoN": "nsp14",
    "NendoU": "nsp15",
    "2'-O-MT": "nsp16",
}

def parse_mut_file(fname):
    ''' Split a list into sub lists, depending on a delimiter. '''
    # 1: read the file
    print(f"Working with: {fname}")
    with open(fname, 'r') as f:
        thelist = [line.rstrip() for line in f.readlines()]
    # 2: get info from filename, set up URLextractor
    bname = os.path.basename(fname)
    gene = bname.split("_")[0]
    # Update gene name if it matches a product
    gene = PRODUCT_TO_NSP.get(gene, gene)  # Convert to nsp if in mapping
    category = '_'.join([x.replace(".txt", "") for x in bname.split("_")[1:]])
    delimiters = ""
    extractor = URLExtract()
    # 3: get info from the file based on: https://stackoverflow.com/questions/45281189/split-list-into-lists-based-on-a-character-occurring-inside-of-an-element
    results = []
    sublist = []
    for item in thelist:
        if item in delimiters:
            results.append(sublist)
            sublist = []
        else:
            sublist.append(item)
    if sublist:
        results.append(sublist)
    results = [x for x in results if len(x) > 0]
    # 4 organise the information
    mut_list = []
    for r in results:
        mut = r.pop()
        flat = ''.join([x.replace("#", "") for x in r])
        extraction = extractor.find_urls(flat)
        reference = ' ; '.join(extraction) if len(extraction) > 0 else 'N/A'
        info = {'gene': gene, 'category': category, 'mutation': mut, 'information': flat, 'reference': reference}
        mut_list.append(info)
    # 5: return the organised dictionary
    return pd.DataFrame(mut_list)


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print(usage)
        sys.exit(1)

    pokay_data_dir = sys.argv[1]
    outfile = sys.argv[2]

    files = glob.glob(pokay_data_dir + "/*.txt")
    df_list = []
    for f in files:
        result = parse_mut_file(f)
        df_list.append(result)
    combined = pd.concat(df_list).drop_duplicates().reset_index(drop=True)
    combined.to_csv(outfile, index=None)
    
    