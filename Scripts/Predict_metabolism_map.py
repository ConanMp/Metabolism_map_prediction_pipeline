#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 24 14:08:41 2020

@author: mconan
"""
import os
import json
import sygma
import argparse
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw

def predict_metabolite(smile, scenario, image_output_dir):
    dict_to_return={"nodes":[],"links":[]}
    parent = Chem.MolFromSmiles(smile)
    metabolic_tree = scenario.run(parent)
    nodes=metabolic_tree.nodes
    source_parent_key=metabolic_tree.parentkey
    source_parent_node=nodes[source_parent_key]
    source_parent_mol=source_parent_node.mol
    source_parent_smile=Chem.MolToSmiles(source_parent_mol)
    image_output=image_output_dir+str(0)+".svg"
    Draw.MolToFile(source_parent_mol, image_output, size=(170, 170), includeAtomNumbers=False)
    Draw.MolToFile(source_parent_mol, image_output.replace(".svg","_atom_num.svg"), size=(170, 170), includeAtomNumbers=True)
    parent_node_dic={"id":0, "Original_smile":smile,"smile_in_sygma_tree":source_parent_smile, "node_key":source_parent_key, "img_svg":"0"}
    dict_to_return["nodes"].append(parent_node_dic)
    id_index=0
    reaction_id=0
    for key_node in nodes.keys():
        if key_node!=source_parent_key:
            id_index+=1
            node=nodes[key_node]
            node_chem_obj=node.mol
            node_smile=Chem.MolToSmiles(node_chem_obj)
            parents=node.parents
            node_image_output=image_output_dir+str(id_index)+".svg"
            Draw.MolToFile(node_chem_obj, node_image_output, size=(170, 170), includeAtomNumbers=False)
            Draw.MolToFile(node_chem_obj, node_image_output.replace(".svg","_atom_num.svg"), size=(170, 170), includeAtomNumbers=True)
            node_dic={"id":id_index,"smile_in_sygma_tree":node_smile, "node_key":key_node, "img_svg":str(id_index)}
            dict_to_return["nodes"].append(node_dic)
            for parent_id in parents:
                reaction_id+=1
                rule_used=parents[parent_id]
                node_parent=nodes[parent_id]
                parent_mol=node_parent.mol
                parent_smile=Chem.MolToSmiles(parent_mol)
                rule_name=rule_used.rulename
                reaction_dic={"id":reaction_id, "source_key":parent_id, "target_key":key_node, "reaction_name":rule_name, "source_smile":parent_smile, "target_smile":node_smile, "atom_number_of_reaction":0}
                for node_json in dict_to_return["nodes"]:
                    if node_json["node_key"]==parent_id:
                        reaction_dic["source"]=node_json["id"]
                    if node_json["node_key"]==key_node:
                        reaction_dic["target"]=node_json["id"]
                dict_to_return["links"].append(reaction_dic)
    return(dict_to_return)

def write_new_json(json_object, output):
    output_file=open(output, 'w')
    output_file.write(json.dumps(json_object, indent=4))
    output_file.close()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("SMILES", help="the SMILES string of the compound to predict metabolism")
    parser.add_argument("Name", help="the name of the compound to predict metabolism, it will be use to name the metabolism map")
    parser.add_argument("Output_dir", help="The path to the output directory whee the metabolisme map will be save", default=False)
    scenario_AHA = sygma.Scenario([
        [sygma.ruleset['phase1'], 2],
        [sygma.ruleset['phase2'], 2]])

    args = parser.parse_args()
    compound_smiles=args.SMILES
    compound_directory=args.Output_dir
    compound_name=args.Name
    img_dir_output=None

    if compound_directory:
        if compound_directory[len(compound_directory)-1]!='/':
            compound_directory=compound_directory+'/'
        if not os.path.exists(compound_directory+compound_name):
            os.makedirs(compound_directory+compound_name+"/Mol_img/")
            img_dir_output=compound_directory+compound_name+"/Mol_img/"
    else:
        if not os.path.exists(compound_name):
            os.makedirs(compound_name+"/Mol_img/")
            img_dir_output=compound_name+"/Mol_img/"

    metabolism_map=predict_metabolite(compound_smiles, scenario_AHA, img_dir_output)
    metabolism_map_output=compound_directory+compound_name+"/"+compound_name+".json"
    write_new_json(metabolism_map, metabolism_map_output)

main()
