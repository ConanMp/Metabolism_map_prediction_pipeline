#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 22 15:16:00 2019

@author: mconan
"""

from selenium import webdriver
from selenium.common.exceptions import NoSuchElementException
import time
import os
import json
import argparse

""" Make analysis W2D """

def make_W2D_results(node, geckodriver):
    a_smile=node["smile_in_sygma_tree"]
    a_smile=a_smile.replace("+","").replace("[O-]","O").replace("O-","O")
    print("SMILES used: "+a_smile)
    url_w2drug="http://www.way2drug.com/SOMP/"
    firefox_webdriver=webdriver.Firefox(executable_path=geckodriver)
    firefox_webdriver.get(url_w2drug)
    firefox_webdriver.find_element_by_id("myHeader3").click()
    firefox_webdriver.find_element_by_id("smi").clear()
    firefox_webdriver.find_element_by_id("smi").send_keys(a_smile)
    firefox_webdriver.find_element_by_css_selector("input[type=\"button\"]").click()

    w2D_results={}
    firefox_webdriver.switch_to.frame("upload_target")
    results_ready=False
    while not results_ready:
        time.sleep(1)
        try:
            firefox_webdriver.find_element_by_xpath("/html/body/div[2]/div[3]/a[1]")
            results_ready=True
        except NoSuchElementException:
            results_ready=False

    for index_for_result_table in range(1,8):
        firefox_webdriver.switch_to.default_content();    # Permet de passer au centre de la page
        firefox_webdriver.switch_to.frame("upload_target")
        link_str="/html/body/div[2]/div[3]/a["+str(index_for_result_table)+"]"
        firefox_webdriver.find_element_by_xpath(link_str).click()
        time.sleep(0.5)
        firefox_webdriver.switch_to.frame("result_data")

        tbody=firefox_webdriver.find_element_by_xpath("/html/body/table/tbody")
        i=0
        enzyme_name=None
        atom_dic={}
        for row in tbody.find_elements_by_xpath("./tr"):
            j=0    
            atom_key=None
            for cells in row.find_elements_by_xpath("./td"):
                string_value=cells.text
                if i==0:
                    enzyme_name=string_value
                elif i==1:
                    continue
                else:
                    if j==0:
                        atom_dic[string_value]={}
                        atom_key=string_value
                    elif j==1:
                        atom_dic[atom_key]["Rank"]=string_value
                    elif j==2:
                        atom_dic[atom_key]["Score"]=string_value
                j+=1
            i+=1
        w2D_results[enzyme_name]=atom_dic
    firefox_webdriver.quit()
    node["way2drug_results"]=w2D_results

""" Make analysis Xenosite """

def make_xenosite_analysis(node_dic, downloadPath, geckodriver, Make_SORs):
    smile_string=node_dic["smile_in_sygma_tree"]
    url_xenosite="https://swami.wustl.edu/xenosite/submit"

    profile = webdriver.FirefoxProfile()
    profile.set_preference('browser.download.folderList', 2)
    profile.set_preference("browser.download.manager.showWhenStarting", False)
    profile.set_preference("browser.download.dir", downloadPath)
    profile.set_preference("browser.helperApps.neverAsk.openFile","text/tab-separated-values, chemical/x-mdl-sdfile, application/zip")
    profile.set_preference("browser.helperApps.neverAsk.saveToDisk","text/tab-separated-values, chemical/x-mdl-sdfile, application/zip")
    profile.set_preference("browser.helperApps.alwaysAsk.force", False)
    profile.set_preference("browser.download.manager.showAlertOnComplete", False)
    profile.set_preference("browser.download.manager.closeWhenDone", False)

    firefox_webdriver=webdriver.Firefox(executable_path=geckodriver, firefox_profile=profile)
    

    option_reactivity="/html/body/div[2]/div[1]/div/div/div/div/form/div/div[1]/table/tbody/tr/td/select/option[7]"
    option_ugt="/html/body/div[2]/div[1]/div/div/div/div/form/div/div[1]/table/tbody/tr/td/select/option[7]"
    option_cyp="/html/body/div[2]/div[1]/div/div/div/div/form/div/div[1]/table/tbody/tr/td/select/option[3]"

    options_dic={"xenositeMetabolismFilesId":option_cyp, "xenositeUgtFilesId":option_ugt}
    if Make_SORs:
        options_dic={"xenositeReactivityFilesId":option_reactivity,"xenositeMetabolismFilesId":option_cyp, "xenositeUgtFilesId":option_ugt}
    for key in options_dic:
        firefox_webdriver.get(url_xenosite)
        option=options_dic[key]
        format_smile="/html/body/div[2]/div[1]/div/div/div/div/form/div/div[3]/table/tbody/tr[4]/td/select/option[2]"
        query_xpath="//*[@id='molecule_input']"
        launch_button="/html/body/div[2]/div[1]/div/div/div/div/form/div/div[4]/input"
        loaded_page_check="/html/body/div[2]/div[1]/div/div/div/div/form/div/div[3]/table/tbody/tr[1]/th"
        loaded=False
        while(not loaded):
            try:
                text_to_check_query_page=firefox_webdriver.find_element_by_xpath(loaded_page_check).text
                if text_to_check_query_page=="Enter Molecules":
                    loaded=True
            except NoSuchElementException:
                time.sleep(0.5)
        firefox_webdriver.find_element_by_xpath(option).click()
        firefox_webdriver.find_element_by_xpath(format_smile).click()
        firefox_webdriver.find_element_by_xpath(query_xpath).clear()
        firefox_webdriver.find_element_by_xpath(query_xpath).send_keys(smile_string)
        firefox_webdriver.find_element_by_xpath(launch_button).click()

        results_availabel=False
        xpath_to_check="/html/body/div[2]/div[1]/div/div[1]/div/div/table/tbody/tr[2]/td"
        while(not results_availabel):
            try:
                text_to_check=firefox_webdriver.find_element_by_xpath(xpath_to_check).text
                if text_to_check=="Results expire 24 hours after creation":
                    results_availabel=True
            except NoSuchElementException:
                try:
                    check_error_server=firefox_webdriver.find_element_by_xpath("/html/body/h1")
                    if check_error_server.text=="Internal Server Error":
                        firefox_webdriver.refresh()
                except NoSuchElementException:
                    time.sleep(3)
                time.sleep(3)
                firefox_webdriver.refresh()

        job_id=firefox_webdriver.current_url
        job_id=job_id.split('/')
        job_id=job_id[5]
        node_dic[key]=job_id
        prediction="/html/body/div[2]/div[1]/div/div[1]/div/div/table/tbody/tr[3]/td/a[1]"
        sdf="/html/body/div[2]/div[1]/div/div[1]/div/div/table/tbody/tr[3]/td/a[2]"
        figure="/html/body/div[2]/div[1]/div/div[1]/div/div/table/tbody/tr[3]/td/a[3]"
        while not os.path.exists(downloadPath+job_id+".molecules.zip"):
            link_download = firefox_webdriver.find_element_by_xpath(figure)
            link_download.click()
            time.sleep(1)
            while os.path.exists(downloadPath+job_id+".molecules.zip.part"):
                time.sleep(0.5)
        while not os.path.exists(downloadPath+job_id+".predictions.tsv"):
            link_download = firefox_webdriver.find_element_by_xpath(prediction)
            link_download.click()
            time.sleep(1)
        while not os.path.exists(downloadPath+job_id):
            link_download = firefox_webdriver.find_element_by_xpath(sdf)
            link_download.click()
            time.sleep(1)
    firefox_webdriver.quit()

def treat_node_reactivity(json_obj, tools_result_dir):
    for node in json_obj["nodes"]:
        node_reactivity=node["xenositeReactivityFilesId"]
        reactivity_file_result=tools_result_dir+node_reactivity+".predictions.tsv"
        result_file=open(reactivity_file_result, 'r')
        result_lines=result_file.readlines()
        first_line=result_lines[0]
        first_line=first_line.replace('\n','').split('\t')
        DNA_index=None
        node_reactivity=0.00
        for elem_index in range(0, len(first_line)):
            elem=first_line[elem_index]
            if elem=="DNA":
                DNA_index=elem_index
        node["dnaAdductSites"]={}
        for line_index in range(2, len(result_lines)):
            line = result_lines[line_index]
            line=line.replace('\n','').split('\t')
            atom_number=line[0]
            atom_number=atom_number.split('.')
            atom_number=atom_number[len(atom_number)-1]
            score_value=float(line[DNA_index])
            if score_value>node_reactivity:
                node_reactivity=score_value
            if score_value>=0.70:
                node["Reactive_to_DNA"]=True
                node["dnaAdductSites"][atom_number]=score_value
        
        if "Reactive_to_DNA" in node.keys():
            reactive_color="Orange"
            node["Reactive_color"]=reactive_color
            for atom in node["dnaAdductSites"].keys():
                score = node["dnaAdductSites"][atom]
                if score >= 0.85:
                    reactive_color="Red"
                    node["Reactive_color"]=reactive_color
                    break
    return(json_obj)

def write_new_json(json_object, output):
    output_file=open(output, 'w')
    output_file.write(json.dumps(json_object, indent=4))
    output_file.close()

def make_analysis_for_json(json_path, geckodriver, Use_Xenosite, Make_SORs):
    json_analysis_output=json_path.replace(".json", "_with_SOM_results.json")
    json_path_splitted=json_path.split('/')
    
    result_dir_path=[]
    for i in range(0, len(json_path_splitted)-1):
        result_dir_path.append(json_path_splitted[i])
    result_dir_path='/'.join(result_dir_path)+'/'
    tools_results=result_dir_path+"toolsResults/"
    if not os.path.exists(tools_results):
        os.mkdir(tools_results)
    
    print("Files from Xenosite will be save at:\n"+str(tools_results))
    json_file=open(json_path, 'r')
    json_str=json_file.read()
    json_obj=json.loads(json_str)
    json_file.close()
    nodes_list=json_obj["nodes"]
    for node in nodes_list:
        node_id=node["id"]
        print("Treating node ID:\t"+str(node_id))
        make_W2D_results(node, geckodriver)
        if Use_Xenosite:
            make_xenosite_analysis(node, tools_results, geckodriver, Make_SORs)
        print("Node ID:\t"+str(node_id)+" treated.")
        write_new_json(json_obj, json_analysis_output)
    
    if Make_SORs and Use_Xenosite:
        Map_treated=treat_node_reactivity(json_obj, tools_results)
        write_new_json(Map_treated, json_analysis_output)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("Metabolism_map", help="Path of the metabolism map to use to get results of SOM predictions")
    parser.add_argument("geckodriver", help="Path of the geckodriver needed to use mozilla by selenium. The geckodriver to use need to be associated to mozilla")
    parser.add_argument("--DNA_reactivity_prediction", help="option to use if you want to get SOR score of metabolites", default=False)
    parser.add_argument("--Use_Xenosite", help="option to use SOM prediction of Xenosite (website often unavailale due to regular maintenance", default=False)
    
    args = parser.parse_args()
    metabolism_map_path=args.Metabolism_map
    geckodriver=args.geckodriver
    Use_Xenosite=args.Use_Xenosite
    get_SORs=args.DNA_reactivity_prediction
    
    make_analysis_for_json(metabolism_map_path, geckodriver, Use_Xenosite, get_SORs)
    
    print("\nAll SOM analysis done.\n")

main()
