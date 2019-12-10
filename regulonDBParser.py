#!/usr/bin/env python
''' Author  : Huy Nguyen
    Program : Parsing operon txt file from regulon db, and write out an operon file that store operon name, and its gene
                the operon will have to have at least k genes, given from user input
    Start   : 07/28/2019
    End     : 07/30/2019
'''

import argparse

def get_arguments():

    parser = argparse.ArgumentParser(description="This program will be parse regulonDB file and create a gene block file that satisfy how many gene an operon has at leask")
    
    parser.add_argument("-i", "--input", default='./E_Coli/regulonDB.txt',
                help="A regulonDB file that store operon info.")
    
    return parser.parse_args()
    
# main function    
if __name__ == "__main__":
    print (1)