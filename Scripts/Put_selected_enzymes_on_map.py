#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 26 17:27:02 2021

@author: mconan
"""
import json
import argparse

def read_selected_enzyme_file(selected_enzyme_file_path):
    general_enzyme_list=["UGTs","NATs","SULTs","GSTs"]
    selected_enzymes_on_reaction_dic={}
    selected_enzyme_file=open(selected_enzyme_file_path, 'r')
    lines = selected_enzyme_file.readlines()
    selected_enzyme_file.close()
    for line_index in range(1,len(lines)):
        line = lines[line_index]
        line=line.replace('\n','').split('\t')
        reaction_id=int(line[0])
        selected_enzyme=line[1]
        if selected_enzyme not in general_enzyme_list:
            selected_enzyme="CYP"+selected_enzyme
        selected_enzymes_on_reaction_dic[reaction_id]=selected_enzyme
    return(selected_enzymes_on_reaction_dic)

def read_json_file(json_path):
    json_file=open(json_path, 'r')
    json_dic=json.loads(json_file.read())
    json_file.close()
    return(json_dic)

def apply_enzyme_on_reactions(map_json_dic, readed_selected_enzyme_dic):
    for reaction in map_json_dic["links"]:
        reaction_id=reaction["id"]
        reaction_selected_enzyme="None"
        if reaction_id in readed_selected_enzyme_dic.keys():
            reaction_selected_enzyme=readed_selected_enzyme_dic[reaction_id]
        #reaction["reaction_name"]=reaction["reaction_name"]+"  /  "+reaction_selected_enzyme
        reaction["selected_enzyme"]=reaction_selected_enzyme

def write_json(json_object, output):
    output_file=open(output, 'w')
    output_file.write(json.dumps(json_object, indent=4))
    output_file.close()
    
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("metabolism_map", help="the path to the file of the .json metabolism map witch have to be annotated with selected enzymes on reactions.")
    parser.add_argument("selected_enzymes", help="the path to the .tsv file which describe for each reaction of the metabolism map which enzyme was selected.")
    parser.add_argument("--output", help="the path to the new map file where the new map with annoted enzymes on reactions.\n\tIf not used the path used will be the path of the map with .json replaced by _with_selected_enzymes.json", default=False)
    args = parser.parse_args()
    selected_enzyme_path=args.selected_enzymes
    metabolism_map_json_path=args.metabolism_map
    output_path=args.output
    if not output_path:
        output_path=metabolism_map_json_path.replace(".json","_with_selected_enzymes.json")
    readed_selected_enzyme_on_prob_calculation=read_selected_enzyme_file(selected_enzyme_path)
    metabolism_map=read_json_file(metabolism_map_json_path)
    apply_enzyme_on_reactions(metabolism_map, readed_selected_enzyme_on_prob_calculation)
    write_json(metabolism_map, output_path)

main()