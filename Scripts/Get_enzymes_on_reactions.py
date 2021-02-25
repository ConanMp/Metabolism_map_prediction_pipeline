#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: mconan
"""

import json
import argparse

def read_phaseI_sygma(phaseI_path):
    phaseI_reaction_name=set([])
    phaseI_smarts_reaction_file=open(phaseI_path, 'r')
    phaseI_lines=phaseI_smarts_reaction_file.readlines()
    for line in phaseI_lines:
        splited_line=line.replace('\n','').split('\t')
        if len(splited_line)==4:
            reaction_name=splited_line[2]
            phaseI_reaction_name.add(reaction_name)
    phaseI_smarts_reaction_file.close()
    return(list(phaseI_reaction_name))

def read_SMIRKS_label_dic(SMIRKS_label_dic_path):
    nats_Accerlys_names=[]
    cyps_Accerlys_names=[]
    ugts_Accerlys_names=[]
    sults_Accerlys_names=[]
    gsts_Accerlys_names=[]

    SMIRKS_label_dic_file=open(SMIRKS_label_dic_path,'r')
    lines=SMIRKS_label_dic_file.readlines()
    for line in lines:
        line=line.replace('\n','')
        splitted_line=line.split('\t')
        rule_name=splitted_line[0]
        enzyme_family=splitted_line[1]
        if enzyme_family=="CYPs":
            cyps_Accerlys_names.append(rule_name)
        elif enzyme_family=="NATs":
            nats_Accerlys_names.append(rule_name)
        elif enzyme_family=="SULTs":
            sults_Accerlys_names.append(rule_name)
        elif enzyme_family=="UGTs":
            ugts_Accerlys_names.append(rule_name)
        elif enzyme_family=="GSTs":
            gsts_Accerlys_names.append(rule_name)
    return(cyps_Accerlys_names,ugts_Accerlys_names,nats_Accerlys_names,gsts_Accerlys_names,sults_Accerlys_names)

def set_enzymes_on_reaction(json_dic, json_path, use_xenosite, CYP_W2D_thresold, CYP_Xenosite_thresold, UGT_W2D_thresold, UGT_Xenosite, SMIRKS_label_dic_path):
    cyps_Accerlys_names, ugts_Accerlys_names, nats_Accerlys_names, gsts_Accerlys_names, sults_Accerlys_names = read_SMIRKS_label_dic(SMIRKS_label_dic_path)
    reactions=json_dic["links"]
    nodes=json_dic["nodes"]

    xenosite_result_path=None
    if use_xenosite:
        json_path_splitted=json_path.split('/')
        xenosite_result_path=json_path_splitted[:len(json_path_splitted)-1]
        xenosite_result_path='/'.join(xenosite_result_path)
        xenosite_result_path=xenosite_result_path+'/toolsResults/'

    for reaction in reactions:
        w2d_enzyme_results={}
        xenosite_enzyme_results={}
        if reaction["reaction_name"] in cyps_Accerlys_names:
            reaction["enzymes"]={}
            atom_number_of_reaction=int(reaction["atom_number_of_reaction"])
            source_id=reaction["source"]
            for node in nodes:
                node_id=node["id"]
                if node_id==source_id:

                    # Treatement of Way2Drug SOM prediction of Cytochrome P450
                    w2d_results=node["way2drug_results"]
                    for enzyme in w2d_results.keys():
                        if enzyme=="UGT" or enzyme=="Standard":
                            continue
                        else:
                            enzyme_result=w2d_results[enzyme]
                            atom_of_reaction_result=enzyme_result[str(atom_number_of_reaction)]
                            w2d_SOM_score=float(atom_of_reaction_result["Score"].replace(',','.'))
                            if w2d_SOM_score >= CYP_W2D_thresold :
                                w2d_enzyme_results[enzyme]=w2d_SOM_score



                    # Treatement of Xenosite SOM prediction of Cytochrome P450
                    if use_xenosite:
                        xenosite_result_file_path=xenosite_result_path+node["xenositeMetabolismFilesId"]+".predictions.tsv"
                        xenosite_cyp_file=open(xenosite_result_file_path,'r')
                        xenosite_cyp_lines=xenosite_cyp_file.readlines()
                        for a_xeno_cyp_line in xenosite_cyp_lines:
                            if a_xeno_cyp_line.startswith("Model"):
                                continue
                            else:
                                a_xeno_cyp_line=a_xeno_cyp_line.replace('\n','')
                                a_xeno_cyp_line=a_xeno_cyp_line.split('\t')
                                model_xeno_cyp_file=a_xeno_cyp_line[0]
                                if model_xeno_cyp_file == "HLM":
                                    continue
                                atom_xeno_cyp_file=int(a_xeno_cyp_line[2])
                                score_xeno_cyp_file=float(a_xeno_cyp_line[5])
                                if atom_xeno_cyp_file==(atom_number_of_reaction+1) and score_xeno_cyp_file>=CYP_Xenosite_thresold:
                                    xenosite_enzyme_results[model_xeno_cyp_file]=score_xeno_cyp_file
                        xenosite_cyp_file.close()
                        break
            if w2d_enzyme_results!={}:
                reaction["enzymes"]["Way2Drug"]=w2d_enzyme_results
            if xenosite_enzyme_results!={}:
                reaction["enzymes"]["Xenosite"]=xenosite_enzyme_results


        elif reaction["reaction_name"] in ugts_Accerlys_names:
            reaction["enzymes"]={}
            atom_number_of_reaction=int(reaction["atom_number_of_reaction"])
            source_id=reaction["source"]
            for node in nodes:
                node_id=node["id"]
                if node_id==source_id:

                    # Treatement of Xenosite SOM prediction of UGT
                    xenosite_result_file_path=xenosite_result_path+node["xenositeUgtFilesId"]+".predictions.tsv"
                    xenosite_ugt_file=open(xenosite_result_file_path,'r')
                    xenosite_ugt_lines=xenosite_ugt_file.readlines()
                    for a_xeno_ugt_line in xenosite_ugt_lines:
                        if a_xeno_ugt_line.startswith("XenoSite"):
                            a_xeno_ugt_line=a_xeno_ugt_line.replace('\n','')
                            a_xeno_ugt_line=a_xeno_ugt_line.split('\t')
                            atom_xeno_ugt_file=int(a_xeno_ugt_line[2])
                            score_xeno_ugt_file=float(a_xeno_ugt_line[3])
                            if atom_xeno_ugt_file==(atom_number_of_reaction+1) and score_xeno_ugt_file>=UGT_Xenosite:
                               xenosite_enzyme_results["UGT"]=score_xeno_ugt_file
                    xenosite_ugt_file.close()

                    # Treatement of Way2Drug SOM prediction of UGT
                    w2d_UGT_results=node["way2drug_results"]["UGT"]
                    w2d_SOM_score=float(w2d_UGT_results[str(atom_number_of_reaction)]["Score"].replace(',','.'))
                    if w2d_SOM_score>=UGT_W2D_thresold:
                        w2d_enzyme_results["UGT"]=w2d_SOM_score
                    break
            if w2d_enzyme_results!={}:
                reaction["enzymes"]["Way2Drug"]=w2d_enzyme_results
            if xenosite_enzyme_results!={}:
                reaction["enzymes"]["Xenosite"]=xenosite_enzyme_results


def write_json(json_dic, output_path):
    output_file=open(output_path, 'w')
    output_file.write(json.dumps(json_dic, indent=4))
    output_file.close()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("json_map_path", help="the .json file path where the metabolic map is saved")
    parser.add_argument("SMIRKS_label_file", help="the path to file where SMIRKS labels are associated to enzymes")
    parser.add_argument("--output", help="the file where you want to write the new json file with enzymes linked to reactions. If you don't use this argument the original file will be overwrite.")
    parser.add_argument("--way2drug_cyp_thresold", default="0.0", help="a thresold for Site Of Metabolism (SOM) score prediction of Cytochrome P450 from Way2Drug If not defined it will bet set at 0.0")
    parser.add_argument("--way2drug_ugt_thresold", default="0.0", help="a thresold for SOM score prediction for UDP-glucuronyl transferase from Way2Drug If not defined it will bet set at 0.00")
    parser.add_argument("--xenosite_cyp_thresold", default="0.0", help="a thresold for SOM score prediction for Cytochrome P450 transferase from Xenosite If not defined it will bet set at 0.00")
    parser.add_argument("--xenosite_ugt_thresold", default="0.0", help="a thresold for SOM score prediction for UDP-glucuronyl transferase from Xenosite If not defined it will bet set at 0.00")
    parser.add_argument("--Use_Xenosite", default=False, help="an argument to use to annotate reactions with xenosite SOM scores")

    args = parser.parse_args()
    SMIRKS_label_dic_path=args.SMIRKS_label_file
    metabolic_map_path=args.json_map_path
    map_output=None
    if args.output==None:
        map_output=args.json_map_path
    else:
        map_output=args.output
    W2D_CYP=float(args.way2drug_cyp_thresold)
    XENO_CYP=float(args.xenosite_cyp_thresold)
    W2D_UGT=float(args.way2drug_ugt_thresold)
    XENO_UGT=float(args.xenosite_ugt_thresold)
    use_xenosite=args.Use_Xenosite

    map_file=open(metabolic_map_path, 'r')    
    map_string=map_file.read()
    map_json=json.loads(map_string)
    map_file.close()
    set_enzymes_on_reaction(map_json, metabolic_map_path, use_xenosite, W2D_CYP, XENO_CYP, W2D_UGT, XENO_UGT, SMIRKS_label_dic_path)
    write_json(map_json, map_output)

main()
