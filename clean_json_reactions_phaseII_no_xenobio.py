#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 23 09:55:10 2019

@author: mconan
"""

import os
import argparse
import time
import itertools
from decimal import Decimal
from ast import literal_eval
import json
import concurrent.futures
from concurrent.futures import ProcessPoolExecutor, as_completed, ThreadPoolExecutor

def get_node_id_list(json_obj):
    node_list=[]
    nodes=json_obj["nodes"]
    for node in nodes:
        node_id=node["id"]
        node_list.append(int(node_id))
    return(node_list)

def get_links_id_list(json_obj):
    node_list=[]
    nodes=json_obj["links"]
    for node in nodes:
        node_id=node["id"]
        node_list.append(int(node_id))
    return(node_list)

def get_reaction_to_del(json_obj):
    reaction_not_supported=["methylation_(aromatic_OH)", "methylation_(thiol)", "glycination_(aromatic_carboxyl)", "glycination_(aliphatic_carboxyl)", "phosphorylation", "dephosphorylation"]
    reaction_id_to_del=[]
    for reaction in json_obj["links"]:
        reaction_name=reaction["reaction_name"]
        if reaction_name in reaction_not_supported:
            reaction_id_to_del.append(reaction["id"])
    return(reaction_id_to_del)

def del_reaction(json_obj, reaction_to_del):
    new_json={"nodes":[], "links":[]}
    new_json["nodes"]=json_obj["nodes"]
    for reaction in json_obj["links"]:
        reaction_id=reaction["id"]
        if reaction_id in reaction_to_del:
            continue
        new_reaction_list=new_json["links"]
        new_reaction_list.append(reaction)
    return(new_json)

def del_N_atom(json_obj):
    new_json={"nodes":[], "links":[]}
    for node_dic in json_obj["nodes"]:
        node_smiles=node_dic["smile_in_sygma_tree"]
        if node_smiles=="N":
            continue
        new_json["nodes"].append(node_dic)
    for reaction in json_obj["links"]:
        source=reaction["source_smile"]
        target=reaction["target_smile"]
        if source == "N":
            continue
        if target == "N":
            continue
        new_json["links"].append(reaction)
    return(new_json)

def del_unrooted_nodes(json_obj):
    new_json={"nodes":[],"links":[]}
    
    rooted_nodes=[0]
    for reaction in json_obj["links"]:
        source=reaction["source"]
        target=reaction["target"]
        if source==0:
            rooted_nodes.append(target)
    for reaction in json_obj["links"]:
        source=reaction["source"]
        target=reaction["target"]
        if source!=0 and source in rooted_nodes:
            rooted_nodes.append(target)

    for node in json_obj["nodes"]:
        node_id=node["id"]
        if node_id not in rooted_nodes:
            continue
        new_nodes_list=new_json["nodes"]
        new_nodes_list.append(node)


    for reaction in json_obj["links"]:
        reaction_source=reaction["source"]
        if reaction_source not in rooted_nodes:
            continue
        new_reactions=new_json["links"]
        new_reactions.append(reaction)

    return(new_json)

def count_numbers(json_obj):
    nodes=json_obj["nodes"]
    links=json_obj["links"]
    return(len(nodes),len(links))

def write_new_json(json_object, output):
    output_file=open(output, 'w')
    output_file.write(json.dumps(json_object, indent=4))
    output_file.close()

Ahas_dir="/home/mconan/Documents/Bayesian_model/AHA_with_sygma/Predict_map_with_sygma/Result/All_AHAs/"

all_haas_path=os.listdir(Ahas_dir)
for AHA_name in all_haas_path:
    Aha_dir=Ahas_dir+AHA_name+"/"
    AHA_path=Aha_dir+AHA_name+"_reduced.json"
    AHA_file=open(AHA_path, 'r')
    AHA_json=json.loads(AHA_file.read())
    AHA_file.close()
    node_id_list=get_node_id_list(AHA_json)
    AHA_json=del_N_atom(AHA_json)
    reactions_to_del=get_reaction_to_del(AHA_json)
    AHA_json=del_reaction(AHA_json, reactions_to_del)
    AHA_json=del_unrooted_nodes(AHA_json)
    num_nodes, num_links=count_numbers(AHA_json)
    print(AHA_name+" : compounds = "+str(num_nodes)+" ; reactions : "+str(num_links))
    AHA_output=AHA_path.replace("_reduced.json","_totaly_reduced.json")
    write_new_json(AHA_json, AHA_output)





#def rename_ids(json_obj):
#    ids_nodes=sorted(get_node_id_list(json_obj))
#    ids_links=sorted(get_links_id_list(json_obj))
#    node_id_equiv_dic={}
#    links_id_equiv_dic={}
#    for i in range (0,len(ids_nodes)):
#        node_id_equiv_dic[ids_nodes[i]]=i
#    for i in range (0, len(ids_links)):
#        links_id_equiv_dic[ids_links[i]]=i
#        
#    for node in json_obj["nodes"]:
#        node_id=node["id"]
#        node["id"]=node_id_equiv_dic[node_id]
#    for reaction in json_obj["links"]:
#        reaction_id=reaction["id"]
#        source=reaction["source"]
#        target=reaction["target"]
#        reaction["id"]=links_id_equiv_dic[reaction_id]
#        reaction["source"]=node_id_equiv_dic[source]
#        reaction["target"]=node_id_equiv_dic[target]
