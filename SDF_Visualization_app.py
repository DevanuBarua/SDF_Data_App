import sys
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Lipinski
from rdkit.Chem import Crippen
from rdkit.Chem import MolSurf
import sqlite3
import tkinter as tk
from tkinter import ttk
from tkinter import filedialog as fd

# Create tkinter main app
window = tk.Tk()
window.title('Bioactive Molecule Data')
window.geometry('1800x900')
window.resizable(False, False)

# Create tkinter windows
frame3 = tk.Frame(window)
frame3.place(height=900,width=1800)

frame1 = tk.LabelFrame(window, text='SDF File')
frame1.place(height=650,width=1700,relx=0.030)

frame2 = tk.LabelFrame(window, text='Functions')
frame2.place(height=200,width=500,rely=0.75,relx=0.35)


# Create buttons for features
Load_Button = ttk.Button(frame2, text="Load Data",command=lambda: load_sdf())
Load_Button.place(height=50,width=120, relx=0.75)

Lipinkski_button = ttk.Button(frame2, text="Lipinski",command=lambda: filter_lipinski())
Lipinkski_button.place(height=50,width=120)

Lead_button = ttk.Button(frame2, text="Lead Likeness", command=lambda: filter_lead())
Lead_button.place(height=50,width=120,rely=0.30)

Bio_button = ttk.Button(frame2, text="Bioavailability", command=lambda: filter_bio())
Bio_button.place(height=50,width=120,rely=0.60)

Save_button = ttk.Button(frame2, text="Save DB", command=lambda: save_db())
Save_button.place(height=50,width=120,relx=0.75,rely=0.30)

Exit_button = ttk.Button(frame2, text="Exit", command=lambda: sys.exit())
Exit_button.place(height=50,width=120,relx=0.75,rely=0.60)

# Creating the base for the data to be shown
columns = ('mol_id','formula','smiles','molwt','Hacceptors','Hdonors','LogP','LogD','ringcount','rotbondcount','tpsa')

mol_data = ttk.Treeview(frame1, columns=columns, show='headings')
mol_data.place(relwidth=1,relheight=1)

# additional widgets for ease of data exploration
yscrollbar = ttk.Scrollbar(frame3,orient=tk.VERTICAL,command=mol_data.yview)
xscrollbar = ttk.Scrollbar(frame3,orient=tk.HORIZONTAL,command=mol_data.xview)
mol_data.configure(yscrollcommand=yscrollbar.set)
mol_data.configure(xscrollcommand=xscrollbar.set)
yscrollbar.place(height=100,width=25,rely=0.3, relx=0.012)
xscrollbar.place(height=20,width=100,rely=0.73,relx=0.46)

# Placing the columns into the app body
mol_data.column('mol_id',anchor='center',width=50)
mol_data.column('formula',anchor='center')
mol_data.column('smiles',anchor='center',width=600)
mol_data.column('molwt',anchor='center',width=100)
mol_data.column('Hacceptors',anchor='center',width=120)
mol_data.column('Hdonors',anchor='center',width=100)
mol_data.column('LogP',anchor='center',width=50)
mol_data.column('LogD',anchor='center',width=50)
mol_data.column('ringcount',anchor='center',width=100)
mol_data.column('rotbondcount',anchor='center',width=120)
mol_data.column('tpsa',anchor='center',width=100)


# Initialization of the columns with their tags
mol_data.heading('mol_id',text="ID")
mol_data.heading('formula',anchor='c',text="Formula" )
mol_data.heading('smiles',text="Smiles" )
mol_data.heading('molwt',text="Mol Weight")
mol_data.heading('Hacceptors',text="H Acceptors")
mol_data.heading('Hdonors',text="H Donors")
mol_data.heading('LogP',text="Log P")
mol_data.heading('LogD',text="Log D")
mol_data.heading('ringcount',text="Ring Count")
mol_data.heading('rotbondcount',text="Rot Bond Count")
mol_data.heading('tpsa',text="TPSA")

# Function for data parsing
def load_sdf():
    file = fd.askopenfilename() # input of sdf file to be opened
    mol_list = []
    with Chem.SDMolSupplier(file) as suppl:
        for mol in suppl:
            each_mol = []
            if mol is not None:
                # Extract relevant data from each relevant record and add them to 'each_mol' list
                molid = mol.GetProp('Mol_ID')
                each_mol.append(molid)

                formula = mol.GetProp('Formula')
                each_mol.append(formula)

                smiles = Chem.MolToSmiles(mol)
                each_mol.append(smiles)

                molwt = Chem.rdMolDescriptors.CalcExactMolWt(mol)
                molwt = round(molwt,2)
                each_mol.append(molwt)

                acceptors = Chem.Lipinski.NumHAcceptors(mol)
                each_mol.append(acceptors)

                donors = Chem.Lipinski.NumHDonors(mol)
                each_mol.append(donors)

                logP = Chem.Crippen.MolLogP(mol)
                logP = round(logP,2)
                each_mol.append(logP)

                logD = mol.GetProp('LogD')
                each_mol.append(logD)

                ringcount = Chem.Lipinski.RingCount(mol)
                each_mol.append(ringcount)

                rotbondcount = Chem.Lipinski.NumRotatableBonds(mol)
                each_mol.append(rotbondcount)

                tpsa = Chem.MolSurf.TPSA(mol)
                tpsa = round(tpsa,2)
                each_mol.append(tpsa)

                # append the 'each_mol' list to 'mol_list'
                # 'mol_list' will have individual lists for each molecule
                mol_list.append(each_mol)

    # add data in 'mol_list' to app
    for row in mol_list:
        mol_data.insert('',tk.END,values=row)


# Function to save data shown on app
def save_db():
    file = fd.asksaveasfilename()
    if not file.endswith(".db"):
        file = file + ".db"
    # Initialize connection to a SQL file
    connection = sqlite3.connect(file)
    c = connection.cursor()
    # Generate empty SQL base
    c.executescript("""
        CREATE Table Molecule
        (MolID Varchar(3),
        Formula Varchar(30),
        Smiles Varchar(30),
        MolWt Int,
        Hacceptors Real,
        Hdonors Real,
        LogP Int,
        LogD Int,
        RingCount Real,
        RotBondCount Real,
        TPSA Int
        )
        """)
    # Insert data from app in SQL database, row by row
    for row_id in mol_data.get_children():
        row = mol_data.item(row_id)['values']
        c.execute("""INSERT INTO Molecule VALUES(?,?,?,?,?,?,?,?,?,?,?)""", row)
    # Save SQL file
    connection.commit()
    connection.close()

# Functions for the different features shown on app
# When pressed, data will be filtered and only qualifying data will be shown
# Remaining data deleted for that instance but will be recovered when process is stopped.
def filter_lipinski():
    filter_list = []
    # Various filters for lipinski
    for item in mol_data.get_children():
        values = mol_data.item(item)['values']
        if float(values[3]) <= 500:
            if float(values[6]) <= 5:
                if int(values[5]) <= 5:
                    if int(values[4]) <= 10:
                        filter_list.append(values)
    # Delete all rows to create empty columns
    mol_data.delete(*mol_data.get_children())
    # repopulate columns with filtered data
    for row2 in filter_list:
        mol_data.insert('','end',text='',values=row2)

# Remaining functions work the same way as Lipinski filter
def filter_lead():
    filter_list = []
    for item in mol_data.get_children():
        values = mol_data.item(item)['values']
        if float(values[3]) <= 450:
            if float(values[7]) <= 4 and float(values[6]) >= -4:
                if int(values[5]) <= 5:
                    if int(values[4]) <= 8:
                        if int(values[8]) <= 4:
                            if int(values[9]) <= 10:
                                filter_list.append(values)
    mol_data.delete(*mol_data.get_children())
    for row2 in filter_list:
        mol_data.insert('','end',text='',values=row2)

def filter_bio():
    filter_list = []
    for item in mol_data.get_children():
        values = mol_data.item(item)['values']
        if float(values[3]) <= 500:
            if float(values[7]) <= 5:
                if int(values[5]) <= 5:
                    if int(values[4]) <= 10:
                        if int(values[9]) <= 10:
                            if float(values[10]) <= 200:
                                filter_list.append(values)
    mol_data.delete(*mol_data.get_children())
    for row2 in filter_list:
        mol_data.insert('', 'end', text='', values=row2)

# initialize app
window.mainloop()