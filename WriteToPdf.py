# -*- coding: utf-8 -*-
#zapiši v pdffajle glavni in pomožni seznam po čredah
import pandas as pd
import numpy as np

glavniSeznam = pd.read_csv('/home/janao/Documents/F4F/OdbiraZivali/RjaveKrave_928_15022017.csv', sep=" ")
glavniSeznam = pd.DataFrame(glavniSeznam, columns=['CRE_SIFRA_CREDA','ID_ZIVALI','DAT_ROJSTVO','ZIV_ID_SEQ'])

pomSeznam = pd.read_csv('/home/janao/Documents/F4F/OdbiraZivali/seznamB_15022017.csv', sep=" ")
pomSeznam = pd.DataFrame(pomSeznam, columns=['CRE_SIFRA_CREDA','ID_ZIVALI','DAT_ROJSTVO','ZIV_ID_SEQ','rel'])

crede = set(glavniSeznam['CRE_SIFRA_CREDA'])

for creda in crede:
    credaDF = glavniSeznam[glavniSeznam['CRE_SIFRA_CREDA'] == creda]
    credaDFb = pomSeznam[pomSeznam['CRE_SIFRA_CREDA'] == creda]
    credaSeznam = pd.concat([credaDF, credaDFb])
    
 
from IPython.display import HTML
h = HTML(d.to_html())
my_file = open('some_file.html', 'w')
my_file.write(h.data)
my_file.close()   

import reportlab
import pyfpdf
from pyfpdf import FPDF, HTMLMixin

class MyFPDF(FPDF, HTMLMixin):
    pass

pdf=MyFPDF()
#First page
pdf.add_page()
pdf.write_html(html)
pdf.output('html.pdf','F')
import pdfkit                
                       
from jinja2 import Environment, FileSystemLoader
env = Environment(loader=FileSystemLoader('.'))
template = env.get_template("myreport.html")