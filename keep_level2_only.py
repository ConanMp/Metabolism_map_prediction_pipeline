#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 20 15:20:44 2020

@author: mconan
"""

import json

import os

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

def rename_ids(json_obj):
    ids_nodes=sorted(get_node_id_list(json_obj))
    ids_links=sorted(get_links_id_list(json_obj))
    node_id_equiv_dic={}
    links_id_equiv_dic={}
    for i in range (0,len(ids_nodes)):
        node_id_equiv_dic[ids_nodes[i]]=i
    for i in range (0, len(ids_links)):
        links_id_equiv_dic[ids_links[i]]=i
        
    for node in json_obj["nodes"]:
        node_id=node["id"]
        node["id"]=node_id_equiv_dic[node_id]
    for reaction in json_obj["links"]:
        reaction_id=reaction["id"]
        source=reaction["source"]
        target=reaction["target"]
        reaction["id"]=links_id_equiv_dic[reaction_id]
        reaction["source"]=node_id_equiv_dic[source]
        reaction["target"]=node_id_equiv_dic[target]

def write_new_json(json_object, output):
    output_file=open(output, 'w')
    output_file.write(json.dumps(json_object, indent=4))
    output_file.close()

All_AHA_dirs="/home/mconan/Documents/Bayesian_model/AHA_with_sygma/Predict_map_with_sygma/Result/All_AHAs/"
all_AHA=os.listdir(All_AHA_dirs)
for AHA in all_AHA:
    AHA_map_path=All_AHA_dirs+AHA+"/"+AHA+".json"
    AHA_file=open(AHA_map_path, 'r')
    AHA_json=json.loads(AHA_file.read())
    AHA_file.close()
    put_source_on_not_sourced(AHA_json)
    reduced_AHA=del_more_than_level2_reactions(AHA_json)
    AHA_map_output=AHA_map_path.replace(".json","_reduced.json")
    write_new_json(reduced_AHA, AHA_map_output)

#PhIP_path="/home/mconan/Documents/Bayesian_model/AHA_with_sygma/Predict_map_with_sygma/Result/PhIP_predicted_map_14_10_20.json"
#MeIQx_path="/home/mconan/Documents/Bayesian_model/AHA_with_sygma/Predict_map_with_sygma/Result/MeIQx_predicted_map_14_10_20.json"

#PhIP_file=open(PhIP_path, 'r')
#PhIP_json=json.loads(PhIP_file.read())
#PhIP_file.close()

#MeIQx_file=open(MeIQx_path, 'r')
#MeIQx_json=json.loads(MeIQx_file.read())
#MeIQx_file.close()

#put_source_on_not_sourced(MeIQx_json)
#put_source_on_not_sourced(PhIP_json)

#reduced_PhIP=del_more_than_level2_reactions(PhIP_json)
#reduced_MeIQx=del_more_than_level2_reactions(MeIQx_json)

#rename_ids(reduced_PhIP)
#rename_ids(reduced_MeIQx)

#MeIQx_outpath="/home/mconan/Documents/Bayesian_model/AHA_with_sygma/Predict_map_with_sygma/Result/MeIQx_predicted_map_reduced_level2_20_10_20.json"
#PhIP_outpath="/home/mconan/Documents/Bayesian_model/AHA_with_sygma/Predict_map_with_sygma/Result/PhIP_predicted_map_reduced_level2_20_10_20.json"

#write_new_json(reduced_MeIQx, MeIQx_outpath)
#write_new_json(reduced_PhIP, PhIP_outpath)
