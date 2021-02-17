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
    
def put_source_on_not_sourced(json_obj):
    new_reactions=[]
    for reaction in json_obj["links"]:
        if "source" not in reaction.keys():
            source_key=reaction["source_key"]
            for node in json_obj["nodes"]:
                node_key=node["node_key"]
                node_id=node["id"]
                if source_key==node_key:
                    reaction["source"]=node_id
        new_reactions.append(reaction)
    json_obj["links"]=new_reactions

def del_more_than_level2_reactions(json_obj):
    new_json_obj={"nodes":[], "links":[]}
    rooted_nodes_level1=set([0])
    for reaction in json_obj["links"]:
        source=reaction["source"]
        if source==0:
            rooted_nodes_level1.add(reaction["target"])
            new_reactions=new_json_obj["links"]
            new_reactions.append(reaction)
            new_json_obj["links"]=new_reactions
    rooted_nodes=set([])
    for reaction in json_obj["links"]:
        source=reaction["source"]
        if source!=0 and source in rooted_nodes_level1:
            rooted_nodes.add(reaction["target"])
            new_reactions=new_json_obj["links"]
            new_reactions.append(reaction)
            new_json_obj["links"]=new_reactions
    rooted_nodes=rooted_nodes.union(rooted_nodes_level1)
    for node in json_obj["nodes"]:
        node_id=node["id"]
        if node_id in rooted_nodes:
            new_nodes=new_json_obj["nodes"]
            new_nodes.append(node)
            new_json_obj["nodes"]=new_nodes
    return(new_json_obj)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("Map", help="the path to the metabolic map to filter")
    parser.add_argument("--Onlylevel2", help="use this argument if you want to keep the metabolic map with only reactions and metabolites of levels 2 or minus", default=False)
    parser.add_argument("--MapCarac", help="use this argument to display caracteristics of metabolic map before and after every reduce", default=False)

    args = parser.parse_args()
    json_path=args.Map
    mol_name=json_path.replace(".json","")
    mol_name=mol_name.split('/')
    mol_name=mol_name[len(mol_name)-1]
    keeplevel2=args.Onlylevel2
    display_carac=args.MapCarac
    
    json_file=open(json_path, 'r')
    metabolic_map=json.loads(json_file.read())    
    json_file.close()    
    put_source_on_not_sourced(metabolic_map)
    if display_carac:
        num_nodes, num_links=count_numbers(metabolic_map)
        print("Before reducing:\n"+mol_name+" : compounds = "+str(num_nodes)+" ; reactions : "+str(num_links))
    only_2levels_metabolic_map=del_more_than_level2_reactions(metabolic_map)
    if display_carac:
        num_nodes, num_links=count_numbers(only_2levels_metabolic_map)
        print("After reducing to only level 2:\n"+mol_name+" : compounds = "+str(num_nodes)+" ; reactions : "+str(num_links))
    
    if keeplevel2:
        level2output=json_path.replace(".json","_only_level2_reduced.json")
        write_new_json(only_2levels_metabolic_map,level2output)
        
    node_id_list=get_node_id_list(only_2levels_metabolic_map)
    metabolic_map_no_azote=del_N_atom(only_2levels_metabolic_map)
    reactions_to_del=get_reaction_to_del(metabolic_map_no_azote)    
    reduced_metabolic_map=del_reaction(metabolic_map_no_azote, reactions_to_del)
    reduced_metabolic_map=del_unrooted_nodes(reduced_metabolic_map)
    
    if display_carac:
        num_nodes, num_links=count_numbers(reduced_metabolic_map)
        print("After reducing levels and reactions associated to xenobiotics metabolism:\n"+mol_name+" : compounds = "+str(num_nodes)+" ; reactions : "+str(num_links))   
    reduced_metabolic_map_output=json_path.replace(".json","_reduced.json")
    write_new_json(reduced_metabolic_map, reduced_metabolic_map_output)
main()
