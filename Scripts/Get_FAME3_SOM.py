#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 18 16:37:39 2019

@author: mconan
"""

# This script add informations of FAME3 SOM prediction score for differents enzymes of phase I and phase II of xenobiotics metabolism.
# The famillies predicted are NATs, SULTs, GSTs, UGTs, and CYPs.

import subprocess
import json
import os
import shutil
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
    
    

def apply_fame_on_json(json_obj, json_path, path_to_fame3, SMIRKS_label_dic_path):
    original_path=os.path.abspath("")
    json_path=os.path.abspath(json_path)
    path_to_fame3=os.path.abspath(path_to_fame3)
    if not os.path.isdir(path_to_fame3):
        path_to_fame3_split=path_to_fame3.split('/')
        path_to_fame3_new=[]
        for elem_index in range(0, len(path_to_fame3_split)-1):
            path_to_fame3_new.append(path_to_fame3_split[elem_index])
        path_to_fame3='/'.join(path_to_fame3_new)
        del path_to_fame3_split
        del path_to_fame3_new
    os.chdir(path_to_fame3)
    fame3_output=get_fame3_output_directory(json_path)
    tools_result_path=fame3_output
    if not os.path.isdir(tools_result_path):
        os.makedirs(tools_result_path)
    
    cyps_Accerlys_names, ugts_Accerlys_names, nats_Accerlys_names, gsts_Accerlys_names, sults_Accerlys_names = read_SMIRKS_label_dic(SMIRKS_label_dic_path)
    
    if not os.path.isdir(tools_result_path+"NATs/"):
        os.mkdir(tools_result_path+"NATs/")
    if not os.path.isdir(tools_result_path+"CYPs/"):
        os.mkdir(tools_result_path+"CYPs/")
    if not os.path.isdir(tools_result_path+"UGTs/"):
        os.mkdir(tools_result_path+"UGTs/")
    if not os.path.isdir(tools_result_path+"GSTs/"):
        os.mkdir(tools_result_path+"GSTs/")
    if not os.path.isdir(tools_result_path+"SULTs/"):
        os.mkdir(tools_result_path+"SULTs/")

    reactions_list=json_obj["links"]
    compounds_list=json_obj["nodes"]

    for reaction in reactions_list:
        reaction_name=reaction["reaction_name"]
        reaction_source=reaction["source"]

        if reaction_name in nats_Accerlys_names:
            for compound in compounds_list:
                if compound["id"]==reaction_source:
                    if "fame3_results" in compound.keys():
                        if "NATs" in compound["fame3_results"]:
                            break
                    compound_smiles=compound["smile_in_sygma_tree"].replace('(','\(').replace(')','\)')
                    output=tools_result_path+"NATs/"
                    bash_cmd="./fame3 -m acetylations -s "+compound_smiles+" -d 5 -o "+output+" -c"
                    bash_cmd=bash_cmd.replace("'","\\'")
                    subprocess.check_output(['bash', '-c', bash_cmd], stderr=subprocess.STDOUT)
                    results_dir=output+"mol_1_1/"
                    if not os.path.isdir(results_dir):
                        print("Error, no directory name \""+results_dir+"\"")
                    if os.path.isdir(results_dir):
                        os.renames(results_dir, output+"compound_"+str(compound["id"]))
                        if os.listdir(output+"compound_"+str(compound["id"])) == []:
                            no_results=open(output+"compound_"+str(compound["id"])+"/mol_1_1_basic.csv", 'w')
                            no_results.write("No probability error from fame3, No results available.")
                            no_results.close()

                    if "fame3_results" in compound.keys():
                        fame_results=compound["fame3_results"]
                        fame_results["NATs"]="compound_"+str(compound["id"])
                    else:
                        compound["fame3_results"]={"NATs":"compound_"+str(compound["id"])}
                    break

        elif reaction_name in cyps_Accerlys_names:  
            for compound in compounds_list:
                if compound["id"]==reaction_source:
                    if "fame3_results" in compound.keys():
                        if "CYPs" in compound["fame3_results"]:
                            break
                    compound_smiles=compound["smile_in_sygma_tree"].replace('(','\(').replace(')','\)')
                    output=tools_result_path+"CYPs/"
                    bash_cmd="./fame3 -m P1 -s "+compound_smiles+" -d 5 -o "+output+" -c"
                    bash_cmd=bash_cmd.replace("'","\\'")
                    subprocess.check_output(['bash', '-c', bash_cmd], stderr=subprocess.STDOUT)
                    results_dir=output+"mol_1_1/"
                    if not os.path.isdir(results_dir):
                        print("Error, no directory name \""+results_dir+"\"")
                    if os.path.isdir(results_dir):
                        os.renames(results_dir, output+"compound_"+str(compound["id"]))
                        if os.listdir(output+"compound_"+str(compound["id"])) == []:
                            no_results=open(output+"compound_"+str(compound["id"])+"/mol_1_1_basic.csv", 'w')
                            no_results.write("No probability error from fame3, No results available.")
                            no_results.close()

                    if "fame3_results" in compound.keys():
                        fame_results=compound["fame3_results"]
                        fame_results["CYPs"]="compound_"+str(compound["id"])
                    else:
                        compound["fame3_results"]={"CYPs":"compound_"+str(compound["id"])}
                    break

        elif reaction_name in ugts_Accerlys_names:
            for compound in compounds_list:
                if compound["id"]==reaction_source:
                    if "fame3_results" in compound.keys():
                        if "UGTs" in compound["fame3_results"]:
                            break
                    compound_smiles=compound["smile_in_sygma_tree"].replace('(','\(').replace(')','\)')
                    output=tools_result_path+"UGTs/"
                    bash_cmd="./fame3 -m glucuronidations -s "+compound_smiles+" -d 5 -o "+output+" -c"
                    bash_cmd=bash_cmd.replace("'","\\'")
                    subprocess.check_output(['bash', '-c', bash_cmd], stderr=subprocess.STDOUT)
                    results_dir=output+"mol_1_1/"
                    if not os.path.isdir(results_dir):
                        print("Error, no directory name \""+results_dir+"\"")
                    if os.path.isdir(results_dir):
                        os.renames(results_dir, output+"compound_"+str(compound["id"]))
                        if os.listdir(output+"compound_"+str(compound["id"])) == []:
                            no_results=open(output+"compound_"+str(compound["id"])+"/mol_1_1_basic.csv", 'w')
                            no_results.write("No probability error from fame3, No results available.")
                            no_results.close()

                    if "fame3_results" in compound.keys():
                        fame_results=compound["fame3_results"]
                        fame_results["UGTs"]="compound_"+str(compound["id"])
                    else:
                        compound["fame3_results"]={"UGTs":"compound_"+str(compound["id"])}
                    break

        elif reaction_name in sults_Accerlys_names:
            for compound in compounds_list:
                if compound["id"]==reaction_source:
                    if "fame3_results" in compound.keys():
                        if "SULTs" in compound["fame3_results"]:
                            break
                    compound_smiles=compound["smile_in_sygma_tree"].replace('(','\(').replace(')','\)')
                    output=tools_result_path+"SULTs/"
                    bash_cmd="./fame3 -m sulfonations -s "+compound_smiles+" -d 5 -o "+output+" -c"
                    bash_cmd=bash_cmd.replace("'","\\'")
                    subprocess.check_output(['bash', '-c', bash_cmd], stderr=subprocess.STDOUT)
                    results_dir=output+"mol_1_1/"
                    if not os.path.isdir(results_dir):
                        print("Error, no directory name \""+results_dir+"\"")
                    if os.path.isdir(results_dir):
                        os.renames(results_dir, output+"compound_"+str(compound["id"]))
                        if os.listdir(output+"compound_"+str(compound["id"])) == []:
                            no_results=open(output+"compound_"+str(compound["id"])+"/mol_1_1_basic.csv", 'w')
                            no_results.write("No probability error from fame3, No results available.")
                            no_results.close()

                    if "fame3_results" in compound.keys():
                        fame_results=compound["fame3_results"]
                        fame_results["SULTs"]="compound_"+str(compound["id"])
                    else:
                        compound["fame3_results"]={"SULTs":"compound_"+str(compound["id"])}
                    break

        elif reaction_name in gsts_Accerlys_names:
            for compound in compounds_list:
                if compound["id"]==reaction_source:
                    if "fame3_results" in compound.keys():
                        if "GSTs" in compound["fame3_results"]:
                            break
                    compound_smiles=compound["smile_in_sygma_tree"].replace('(','\(').replace(')','\)')
                    output=tools_result_path+"GSTs/"
                    bash_cmd="./fame3 -m GSH-conjugations -s "+compound_smiles+" -d 5 -o "+output+" -c"
                    bash_cmd=bash_cmd.replace("'","\\'")
                    subprocess.check_output(['bash', '-c', bash_cmd], stderr=subprocess.STDOUT)
                    results_dir=output+"mol_1_1/"
                    if not os.path.isdir(results_dir):
                        print("Error, no directory name \""+results_dir+"\"")
                    if os.path.isdir(results_dir):
                        os.renames(results_dir, output+"compound_"+str(compound["id"]))
                        if os.listdir(output+"compound_"+str(compound["id"])) == []:
                            no_results=open(output+"compound_"+str(compound["id"])+"/mol_1_1_basic.csv", 'w')
                            no_results.write("No probability error from fame3, No results available.")
                            no_results.close()

                    if "fame3_results" in compound.keys():
                        fame_results=compound["fame3_results"]
                        fame_results["GSTs"]="compound_"+str(compound["id"])
                    else:
                        compound["fame3_results"]={"GSTs":"compound_"+str(compound["id"])}
                    break
    os.chdir(original_path)

def write_json(json_dic, output_path):
    output_file=open(output_path, 'w')
    output_file.write(json.dumps(json_dic, indent=4))
    output_file.close()

def treat_fame3_results(json_dic, json_path, SMIRKS_label_dic_path):
    tools_result_path=get_fame3_output_directory(json_path)
    reactions_list=json_dic["links"]
    compounds_list=json_dic["nodes"]
    
    cyps_Accerlys_names, ugts_Accerlys_names, nats_Accerlys_names, gsts_Accerlys_names, sults_Accerlys_names = read_SMIRKS_label_dic(SMIRKS_label_dic_path)

    for reaction in reactions_list:
        reaction_name=reaction["reaction_name"]
        reaction_source=reaction["source"]

        if reaction_name in nats_Accerlys_names:
            for compound in compounds_list:
                if compound["id"]==reaction_source:
                    compound_fame_results=compound["fame3_results"]
                    if "NATs" in compound_fame_results.keys():
                        fame3_compound_results_path=tools_result_path+"NATs/"+compound_fame_results["NATs"]+"/mol_1_1_basic.csv"
                        fame3_compound_results_path=fame3_compound_results_path.replace("'","\'")
                        fame3_result_file=open(fame3_compound_results_path, 'r')
                        lines=fame3_result_file.readlines()
                        if len(lines)<=2:
                            if "enzymes" in reaction.keys():
                                reaction["enzymes"]["FAME3"]=0.0
                            else:
                                reaction["enzymes"]={"FAME3":0.0}
                        else:
                            lines=lines[1:len(lines)]
                            atom_number_of_reaction=int(reaction["atom_number_of_reaction"])
                            for line_index in range(0, len(lines)):
                                line=lines[line_index]
                                line=line.split(',')
                                score=float(line[3])
                                atom=str(line[1])
                                atom=atom.split('.')
                                atom=atom[len(atom)-1]
                                atom=int(atom)
                                if atom==(atom_number_of_reaction):
                                    if "enzymes" in reaction.keys():
                                        reaction["enzymes"]["FAME3"]=score
                                    else:
                                        reaction["enzymes"]={"FAME3":score}
                        fame3_result_file.close()

        if reaction_name in cyps_Accerlys_names:
            for compound in compounds_list:
                if compound["id"]==reaction_source:
                    compound_fame_results=compound["fame3_results"]
                    if "CYPs" in compound_fame_results.keys():
                        fame3_compound_results_path=tools_result_path+"CYPs/"+compound_fame_results["CYPs"]+"/mol_1_1_basic.csv"
                        fame3_compound_results_path=fame3_compound_results_path.replace("'","\'")
                        fame3_result_file=open(fame3_compound_results_path, 'r')
                        lines=fame3_result_file.readlines()
                        if len(lines)<=2:
                            if "enzymes" in reaction.keys():
                                reaction["enzymes"]["FAME3"]=0.0
                            else:
                                reaction["enzymes"]={"FAME3":0.0}
                        else:
                            lines=lines[1:len(lines)]
                            atom_number_of_reaction=int(reaction["atom_number_of_reaction"])
                            for line_index in range(0, len(lines)):
                                line=lines[line_index]
                                line=line.split(',')
                                score=float(line[3])
                                atom=str(line[1])
                                atom=atom.split('.')
                                atom=atom[len(atom)-1]
                                atom=int(atom)
                                if atom==(atom_number_of_reaction):
                                    if "enzymes" in reaction.keys():
                                        reaction["enzymes"]["FAME3"]=score
                                    else:
                                        reaction["enzymes"]={"FAME3":score}
                        fame3_result_file.close()

        if reaction_name in sults_Accerlys_names:
            for compound in compounds_list:
                if compound["id"]==reaction_source:
                    compound_fame_results=compound["fame3_results"]
                    if "SULTs" in compound_fame_results.keys():
                        fame3_compound_results_path=tools_result_path+"SULTs/"+compound_fame_results["SULTs"]+"/mol_1_1_basic.csv"
                        fame3_compound_results_path=fame3_compound_results_path.replace("'","\'")
                        fame3_result_file=open(fame3_compound_results_path, 'r')
                        lines=fame3_result_file.readlines()
                        if len(lines)<=2:
                            if "enzymes" in reaction.keys():
                                reaction["enzymes"]["FAME3"]=0.0
                            else:
                                reaction["enzymes"]={"FAME3":0.0}
                        else:
                            lines=lines[1:len(lines)]
                            atom_number_of_reaction=int(reaction["atom_number_of_reaction"])
                            for line_index in range(0, len(lines)):
                                line=lines[line_index]
                                line=line.split(',')
                                score=float(line[3])
                                atom=str(line[1])
                                atom=atom.split('.')
                                atom=atom[len(atom)-1]
                                atom=int(atom)
                                if atom==(atom_number_of_reaction):
                                    if "enzymes" in reaction.keys():
                                        reaction["enzymes"]["FAME3"]=score
                                    else:
                                        reaction["enzymes"]={"FAME3":score}
                        fame3_result_file.close()

        if reaction_name in ugts_Accerlys_names:
            for compound in compounds_list:
                if compound["id"]==reaction_source:
                    compound_fame_results=compound["fame3_results"]
                    if "UGTs" in compound_fame_results.keys():
                        fame3_compound_results_path=tools_result_path+"UGTs/"+compound_fame_results["UGTs"]+"/mol_1_1_basic.csv"
                        fame3_compound_results_path=fame3_compound_results_path.replace("'","\'")
                        fame3_result_file=open(fame3_compound_results_path, 'r')
                        lines=fame3_result_file.readlines()
                        if len(lines)<=2:
                            if "enzymes" in reaction.keys():
                                reaction["enzymes"]["FAME3"]=0.0
                            else:
                                reaction["enzymes"]={"FAME3":0.0}
                        else:
                            lines=lines[1:len(lines)]
                            atom_number_of_reaction=int(reaction["atom_number_of_reaction"])
                            for line_index in range(0, len(lines)):
                                line=lines[line_index]
                                line=line.split(',')
                                score=float(line[3])
                                atom=str(line[1])
                                atom=atom.split('.')
                                atom=atom[len(atom)-1]
                                atom=int(atom)
                                if atom==(atom_number_of_reaction):
                                    if "enzymes" in reaction.keys():
                                        reaction["enzymes"]["FAME3"]=score
                                    else:
                                        reaction["enzymes"]={"FAME3":score}
                        fame3_result_file.close()

        if reaction_name in gsts_Accerlys_names:
            for compound in compounds_list:
                if compound["id"]==reaction_source:
                    compound_fame_results=compound["fame3_results"]
                    if "GSTs" in compound_fame_results.keys():
                        fame3_compound_results_path=tools_result_path+"GSTs/"+compound_fame_results["GSTs"]+"/mol_1_1_basic.csv"
                        fame3_compound_results_path=fame3_compound_results_path.replace("'","\'")
                        fame3_result_file=open(fame3_compound_results_path, 'r')
                        lines=fame3_result_file.readlines()
                        if len(lines)<=2:
                            if "enzymes" in reaction.keys():
                                reaction["enzymes"]["FAME3"]=0.0
                            else:
                                reaction["enzymes"]={"FAME3":0.0}
                        else:
                            lines=lines[1:len(lines)]
                            atom_number_of_reaction=int(reaction["atom_number_of_reaction"])
                            for line_index in range(0, len(lines)):
                                line=lines[line_index]
                                line=line.split(',')
                                score=float(line[3])
                                atom=str(line[1])
                                atom=atom.split('.')
                                atom=atom[len(atom)-1]
                                atom=int(atom)
                                if atom==(atom_number_of_reaction):
                                    if "enzymes" in reaction.keys():
                                        reaction["enzymes"]["FAME3"]=score
                                    else:
                                        reaction["enzymes"]={"FAME3":score}
                        fame3_result_file.close()

def get_map_directory(json_path):
    haa_directory=None
    json_path_splited=json_path.split('/')
    haa_directory=json_path_splited[:len(json_path_splited)-1]
    haa_directory='/'.join(haa_directory)
    haa_directory=haa_directory+'/'
    return(haa_directory)

def get_mol_name(json_path):
    json_path_splited=json_path.split('/')
    mol_name=json_path_splited[len(json_path_splited)-1]
    mol_name=mol_name.split("_")
    mol_name=mol_name[0]
    return(mol_name)
    

def get_fame3_output_directory(HAA_map_path):
    HAA_splited=HAA_map_path.split('/')
    HAA_directory=get_map_directory(HAA_map_path)
    HAA_tools_results=HAA_directory+"toolsResults/"
    if not os.path.exists(HAA_tools_results):
        os.makedirs(HAA_tools_results)
    Haa_name=get_mol_name(HAA_map_path)
    fame3_output_dir="fame3"+HAA_splited[len(HAA_splited)-1].replace('.json','').replace(Haa_name,"")
    fame3_output_dir=HAA_tools_results+fame3_output_dir
    fame3_output_dir=fame3_output_dir+"/"
    return fame3_output_dir


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("map_file", help="the .json file path where the metabolic map is saved")
    parser.add_argument("Fame3", help="the path to the fame3 executable")
    parser.add_argument("SMIRKS_dictionnary", help="the path to the file of the SMIRKS rules labels of SyGMa associated to enzyme familly")
    args = parser.parse_args()
    path_to_json_file=args.map_file
    path_to_directory=get_map_directory(path_to_json_file)
    path_to_fame3_predictor=args.Fame3
    SMIRKS_dictionnary=args.SMIRKS_dictionnary
    
    
    map_json_file=open(path_to_json_file, 'r')
    map_string=map_json_file.read()
    map_json_dic=json.loads(map_string)
    answer=None
    if (path_to_json_file.replace('.json','_fame3_processed.json') in path_to_json_file):
        answer=input("A file already processed by fame3 was found in your directory, do you want to overwrite the file and fame3 results? YES/NO\n\t")
        while (answer!="YES" and answer!="NO"):
            answer=input("A file already processed by fame3 was found in your directory, do you want to overwrite the file and fame3 results? YES/NO\n\t")
    else:
        for file_name in os.listdir(path_to_haa_directory):
            if "_fame3_processed.json" in file_name:
                answer=input("A file already processed by fame3 was found in your directory, do you want to overwrite the file and fame3 results? YES/NO\n\t")
                while (answer!="YES" and answer!="NO"):
                    answer=input("A file already processed by fame3 was found in your directory, do you want to overwrite the file and fame3 results? YES/NO\n\t")

    if answer==None:
        apply_fame_on_json(map_json_dic, path_to_json_file, path_to_fame3_predictor, SMIRKS_dictionnary)
        treat_fame3_results(map_json_dic, path_to_json_file, SMIRKS_dictionnary)
        write_json(map_json_dic, path_to_json_file.replace('.json','_fame3_processed.json'))

    else:
        if answer=="YES":
            fame3_tools=get_fame3_output_directory(path_to_json_file)
            if os.path.isdir(fame3_tools):
                shutil.rmtree(fame3_tools)
                del fame3_tools
            apply_fame_on_json(map_json_dic, path_to_json_file, path_to_fame3_predictor, SMIRKS_dictionnary)
            treat_fame3_results(map_json_dic, path_to_json_file, SMIRKS_dictionnary)
            write_json(map_json_dic, path_to_json_file.replace('.json','_fame3_processed.json'))
    map_json_file.close()

main()
