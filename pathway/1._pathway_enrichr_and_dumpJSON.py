# %pip install --user --upgrade git+https://github.com/Maayanlab/maayanlab-bioinformatics.git

# from maayanlab_bioinformatics.enrichment import enrich_crisp

import pandas as pd
import json
import requests
import time
import glob, os


significance_value = 0.05

preceeding_path = "diff_expr/output"
output_path = "pathway/output"
json_output_path = "pathway/visual"

# Error handling
class NoResults(Exception):
    pass 
class APIFailure(Exception):
    pass

def Enrichr_API(enrichr_gene_list, all_libraries):

    results_df = pd.DataFrame()

    for library_name in all_libraries : 
        ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/addList'
        genes_str = '\n'.join(enrichr_gene_list)
        description = ''
        payload = {
            'list': (None, genes_str),
            'description': (None, description)
        }

        response = requests.post(ENRICHR_URL, files=payload)
        if not response.ok:
            raise APIFailure

        data = json.loads(response.text)
        time.sleep(0.5)
        ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/enrich'
        response = requests.get(f"{ENRICHR_URL}?userListId={data['userListId']}&backgroundType={library_name}")
        if not response.ok:
            raise APIFailure

        data = json.loads(response.text)

        if len(data[library_name]) == 0:
            raise NoResults
        
        res_df  = pd.DataFrame(data[library_name])
        # adds library name to the data frame so the libraries can be distinguished
        res_df['library'] = library_name.replace('_', '')
        results_df = pd.concat([results_df, res_df], axis=0)

    return results_df

def get_batch_enrich(genes, all_libraries, file_name):
    results = Enrichr_API(genes, all_libraries)
    if len(results) > 0:
        res_df = results.rename(columns={
            0: 'rank',
            1: 'term',
            2: 'p-value',
            3: 'zscore',
            4: 'combined_score',
            5: 'overlap_genes',
            6: 'q-value'
        })
        sorted_res_df = res_df.sort_values(by=['library','q-value'], ascending=True)
        filtered_res_df = sorted_res_df[sorted_res_df['q-value'] <= significance_value].reset_index()
        
        filtered_res_df.to_csv(f"{output_path}/{file_name}_enrichment.csv", index=False)
            

def dump_to_json(library_title, compare_group):
    all_tables = []
    batch_pairs = []
    for fullpath_file in glob.glob((f"{output_path}/{compare_group}_*_enrichment.csv")):
        table1 = pd.read_csv(fullpath_file)
        if table1.shape[0]:
            all_tables.append(table1)
            filename = os.path.basename(fullpath_file)
            filename = filename.replace(f"{compare_group}_", "").replace("_enrichment.csv", "")
            batch_pairs.append(filename)

    # Process data to create nodes and links for the graph
    nodes = []
    links = []
    for i, table in enumerate(all_tables):
        for idx, row in enumerate(table.itertuples(index=False)):
            pathway, genes = row.term, eval(row.overlap_genes)
            nodes.append({"id": f"{pathway}({batch_pairs[i]})", 
                          "value": len(genes), 
                          "group": batch_pairs[i]})
            # within table
            for other_idx, other_row in enumerate(table.iloc[idx+1:,].itertuples(index=False)):
                other_pathway, other_genes = other_row.term, eval(other_row.overlap_genes)
                if pathway != other_pathway:
                    shared_genes = set(genes).intersection(set(other_genes))
                    links.append({"source": f"{pathway}({batch_pairs[i]})",
                                "target": f"{other_pathway}({batch_pairs[i]})",
                                "value": len(shared_genes)})
            # cross table
            for j, other_table in enumerate(all_tables[i+1:]):
                for other_idx, other_row in enumerate(other_table.itertuples(index=False)):
                    other_pathway, other_genes = other_row.term, eval(other_row.overlap_genes)
                    shared_genes = set(genes).intersection(set(other_genes))
                    links.append({"source": f"{pathway}({batch_pairs[i]})",
                                "target": f"{other_pathway}({batch_pairs[i+j+1]})",
                                "value": len(shared_genes)})

    with open(f'{json_output_path}/graph_data_{compare_group}.json', 'w') as f:
        # write the data to the file in pretty format
        json.dump({"nodes": nodes,
                "links": links, 
                "groups": [{"name": i} for i in batch_pairs], 
                "title": library_title}, f, indent=4)


# gene_bank_library = ['HDSigDB_Mouse_2021', 'KEGG_2019_Mouse', 'KOMP2_Mouse_Phenotypes_2022', 'Mouse_Gene_Atlas', 'WikiPathways_2019_Mouse']
gene_bank_library = 'KEGG_2019_Mouse'
# for filename in glob.glob(f"{preceeding_path}/*.csv"):
#     gene_list_df = pd.read_csv(filename)
#     genes = gene_list_df.gene.tolist()
#     file_name = os.path.splitext(os.path.basename(filename))[0]
#     if not os.path.exists(f"{output_path}/{file_name}_enrichment.csv"):
#         if len(genes):
#             get_batch_enrich(genes, [gene_bank_library], file_name)


dump_to_json(gene_bank_library, "doseDiff_longTerm")
dump_to_json(gene_bank_library, "doseDiff_shortTerm")
dump_to_json(gene_bank_library, "doseRemoval_expressReturn")
dump_to_json(gene_bank_library, "doseRemoval_longTerm")
dump_to_json(gene_bank_library, "termDiff_doseSame")