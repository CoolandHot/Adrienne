# %pip install --user --upgrade git+https://github.com/Maayanlab/maayanlab-bioinformatics.git

# from maayanlab_bioinformatics.enrichment import enrich_crisp

from typing import Tuple
import pandas as pd
import json
import requests
import time
import glob, os


significance_value = 0.05

preceeding_path = "diff_expr/output"
output_path = "pathway/output"
json_output_path = "visualisation"
colors = ['#5470c6', '#91cc75', '#fac858', '#ee6666', '#73c0de', '#3ba272', '#fc8452', '#9a60b4', '#ea7ccc']

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
            

def node_link_merge_pathways(all_tables: list, batch_pairs: list) -> Tuple[list, list]:
    name_color = {k: v for k, v in zip(batch_pairs, colors[:len(batch_pairs)])}

    # merge pathways with the same name across all tables
    merged_df = pd.concat(all_tables).groupby('term').agg({'overlap_genes': 'sum', 'source': list, 'source_geneLen': list})
    merged_df['overlap_genes'] = merged_df['overlap_genes'].apply(lambda x: list(set(x)))
    merged_df['gene_len'] = merged_df['overlap_genes'].apply(lambda x: len(x))
    merged_df['proportion'] = merged_df[['source_geneLen', 'overlap_genes']].apply(lambda ls: [round(x/len(ls[1]),4) for x in ls[0]], axis=1)
    merged_df = merged_df.drop(['source_geneLen'], axis=1)
    merged_df.reset_index(inplace=True)

    nodes = []
    for index, row in merged_df.iterrows():
        item =  {
            "name": row['term'],
            "value": len(row['overlap_genes']),
            'category': 'mix',
            "itemStyle": {
                "color": {
                    "type": "linear",
                    "colorStops": []
                }
            }
        }
        colorStops = []
        offset = 0
        proportion, source = row['proportion'], row['source']
        print_props = []
        for src, prop in zip(source, proportion):
            colorStops.append({'offset': offset, 'color': name_color[src]})
            p = round(prop/sum(proportion), 3)
            offset += p
            if offset > 1.0:
                offset = 1.0
            colorStops.append({'offset': offset, 'color': name_color[src]})
            print_props.append(f"{src}: {p:.2%}")

        if len(proportion) == 1:
            item['category'] = source[0]
                
        item['itemStyle']['color']['colorStops'] = colorStops
        item['proportion'] = "<br>".join(print_props)
        if item['proportion']=="":
            item['proportion'] = f"{source[0]} 100%"
        nodes.append(item)

    links = []
    for row_idx, row in merged_df.iterrows():
        for other_idx, other_row in merged_df.iloc[row_idx+1:,].iterrows():
            if row['term'] != other_row['term']:
                shared_genes = set(row['overlap_genes']).intersection(set(other_row['overlap_genes']))
                if len(shared_genes)>0:
                    links.append({"source": f"{row['term']}",
                                "target": f"{other_row['term']}",
                                "value": len(shared_genes)})
    return nodes, links

def node_link_separate_sets(all_tables: list, batch_pairs: list) -> Tuple[list, list]:
    nodes = []
    links = []
    for table_idx, table in enumerate(all_tables):
        for row_idx, row in enumerate(table.itertuples(index=False)):
            pathway, genes = row.term, row.overlap_genes
            nodes.append({"name": f"{pathway}({batch_pairs[table_idx]})", 
                          "value": len(genes), 
                          "category": batch_pairs[table_idx]})
            # within table
            for other_idx, other_row in enumerate(table.iloc[row_idx+1:,].itertuples(index=False)):
                other_pathway, other_genes = other_row.term, other_row.overlap_genes
                if pathway != other_pathway:
                    shared_genes = set(genes).intersection(set(other_genes))
                    if len(shared_genes) > 0:
                        links.append({"source": f"{pathway}({batch_pairs[table_idx]})",
                                    "target": f"{other_pathway}({batch_pairs[table_idx]})",
                                    "value": len(shared_genes)})
            # cross table
            for j, other_table in enumerate(all_tables[table_idx+1:]):
                for other_idx, other_row in enumerate(other_table.itertuples(index=False)):
                    other_pathway, other_genes = other_row.term, other_row.overlap_genes
                    shared_genes = set(genes).intersection(set(other_genes))
                    if len(shared_genes) > 0:
                        links.append({"source": f"{pathway}({batch_pairs[table_idx]})",
                                    "target": f"{other_pathway}({batch_pairs[table_idx+j+1]})",
                                    "value": len(shared_genes)})
    return nodes, links

def dump_to_json(library_title, purpose_group):
    all_tables = []
    sample_compare_groups = []
    for fullpath_file in glob.glob((f"{output_path}/{purpose_group}_*_enrichment.csv")):
        table1 = pd.read_csv(fullpath_file)
        if table1.shape[0]:
            filename = os.path.basename(fullpath_file)
            filename = filename.replace(f"{purpose_group}_", "").replace("_enrichment.csv", "")
            sample_compare_groups.append(filename)

            table1['source'] = filename
            table1.overlap_genes = table1.overlap_genes.apply(lambda x: eval(x))
            table1['source_geneLen'] = table1.overlap_genes.apply(lambda x: len(x))
            all_tables.append(table1)

    if len(all_tables) >0:
        # Process data to create nodes and links for the graph
        nodes, links = node_link_separate_sets(all_tables, sample_compare_groups)
        with open(f'{json_output_path}/{purpose_group}_separateSets.json', 'w') as f:
            # write the data to the file in pretty format
            json.dump({"nodes": nodes,
                    "links": links, 
                    "categories": [{"name": i} for i in sample_compare_groups], 
                    "title": library_title}, f, indent=4)
            
        nodes, links = node_link_merge_pathways(all_tables, sample_compare_groups)
        with open(f'{json_output_path}/{purpose_group}_mergePathway.json', 'w') as f:
            # write the data to the file in pretty format
            json.dump({"nodes": nodes,
                    "links": links, 
                    "categories": [{"name": i} for i in sample_compare_groups] + [{'name':"mix","itemStyle":{"color":{"type":"linear","colorStops":[{"offset":0,"color":"#5470c6"},
                                                                                                                                        {"offset":0.3,"color":"#5470c6"},
                                                                                                                                        {"offset":0.3,"color":"#91cc75"},
                                                                                                                                        {"offset":0.6,"color":"#91cc75"},
                                                                                                                                        {"offset":0.6,"color":"#fac858"},
                                                                                                                                        {"offset":1.0,"color":"#fac858"}]}},}], 
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