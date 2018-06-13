import os
from pandas import DataFrame
import numpy as np
import rdkit
from rdkit import DataStructs
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import Draw
from rdkit.Chem import rdFMCS
from rdkit.Chem import AllChem
from rdkit.Chem.Fingerprints import FingerprintMols

TestMol = Chem.MolFromSmiles('CC1OC(=NC1C1=NC(CS1)C1=NC(CS1)C(=O)O)c1ccccc1O')
TestMols=[("Example", TestMol)]

for TestMolName, TestMol in TestMols:
     print "Starting Mol: ", TestMolName
     
     MyMol=Chem.RemoveHs(TestMol)
     R=MyMol.GetRingInfo()
     NumRings = R.NumRings()
     
     # useHs=False,branchedPaths=False
     MyMolFp = Chem.RDKFingerprint(MyMol,fpSize=2048, minPath=7, useBondOrder=False, nBitsPerHash=1, minSize=2048, useHs=False)
     MyMolList = [MyMol]
     #####################################################################################################
     ###############################  Start Query Structure Preparation  #####################
     
     nm=[]
     print " Number of Rings is: ", NumRings
     if( NumRings > 6): NumRings = 6
     print " Number of loops is: ", NumRings
     for rings in range(NumRings):
         AList =[]
         for AMol in MyMolList:
             nm = []
             Chem.MolToSmiles(AMol,kekuleSmiles=True)
             AMol.UpdatePropertyCache(strict=False)
             AllChem.Cleanup(AMol)
             AllChem.SetAromaticity(AMol)
             if AMol.HasSubstructMatch(Chem.MolFromSmarts('[r5,r6,R5,R6]@[r5,r6,R5,R6]')):
                 bis = AMol.GetSubstructMatches(Chem.MolFromSmarts('[r5,r6,R5,R6]@[r5,r6,R5,R6]'))
                 if (len(bis) > 12): bis = bis[0:2]
                 for Beg, End in bis:
                     b = [AMol.GetBondBetweenAtoms(Beg,End).GetIdx()]
                     m =  Chem.FragmentOnBonds(AMol,b, addDummies=False)
                     m.UpdatePropertyCache(strict=False)
                     AllChem.Cleanup(m)
                     AllChem.SetAromaticity(m)    
                     Chem.SanitizeMol(m, (Chem.SanitizeFlags.SANITIZE_ADJUSTHS | Chem.SanitizeFlags.SANITIZE_SYMMRINGS))
                     Chem.SanitizeMol(m, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL^Chem.SanitizeFlags.SANITIZE_KEKULIZE)
                     m.UpdatePropertyCache(strict=False)
                     Chem.FastFindRings(m)
                     AllChem.Cleanup(m)
                     AllChem.SetAromaticity(m)
                     Chem.SetHybridization(m)
                     AList.append(m)
         MyMolList.extend(AList)
     
     uniqps = {}
     for p in MyMolList:
        smi = Chem.MolToSmiles(p)
        uniqps[smi] = p
     MyMolList = uniqps.values()
     
     print len(MyMolList)
     print Chem.MolToSmiles(MyMol)
     #####################################################
     
     rxns=[]
     rxns.append(AllChem.ReactionFromSmarts("[#6:1]-[#8:2]-[CH2:3]>>[#6:1]-[#8H:2].[C:3]"))
     rxns.append(AllChem.ReactionFromSmarts("[#6:1]-[N;!R:2](-[CH3:3])(-[CH3:4])>>([#6:1]-[NH2:2]).([CH4:3]).([CH4:4])"))
     rxns.append(AllChem.ReactionFromSmarts("[#6:1]-[#7:2](-[CH3:4])-[#6:3]>>([#6:1]-[#7H:2]-[#6:3]).([CH4:4])"))
     
     # remove Hydroxy
     rxns.append(AllChem.ReactionFromSmarts("[#6:1]-[OX2H]>>[#6:1].[OH2]"))
     # remove epoxides
     rxns.append(AllChem.ReactionFromSmarts("[C;r3:1][O;r3][C;r3:2]>>[C:1]=[C:2].[OH2]"))
     # remove methoxy
     rxns.append(AllChem.ReactionFromSmarts("[C:1][O][CH3:2]>>[C:1].[CH3OH]"))
     
     
     pattD = Chem.MolFromSmarts('[#6]-[#8]')
     pattR = Chem.MolFromSmarts('[#6]=[#8]')
     replR = Chem.MolFromSmarts('[#6]-[#8]')
     
     print " The number of molecules is now: ", len(MyMolList)
     
     mols = []
     for AMol in MyMolList:
         mL=[]
     
         for rxn in rxns:
             mr = rxn.RunReactants((AMol, ))
             ms = set([Chem.MolToSmiles(x[0]) for x in mr])
             for Rot,Mol in enumerate(ms):
                m = Chem.MolFromSmiles(Mol,sanitize=False)
                AllChem.Cleanup(m)
                Chem.SanitizeMol(m, (Chem.SanitizeFlags.SANITIZE_ADJUSTHS | Chem.SanitizeFlags.SANITIZE_SYMMRINGS))
                if Descriptors.MolWt(m) > 150: mL.append(m)
             mols.extend(mL)
     
         m = AllChem.DeleteSubstructs(AMol,pattD)
         Chem.MolToSmiles(m,kekuleSmiles=True)
         AllChem.Cleanup(m)
         AllChem.SetAromaticity(m)
         mols.append(m)
     
         mL = []
         m = AllChem.ReplaceSubstructs(AMol,pattR,replR)
         for Rot,Mol in enumerate(m):
             Mol.UpdatePropertyCache(strict=False)
             AllChem.Cleanup(Mol)
             AllChem.SetAromaticity(Mol)
             if Descriptors.MolWt(Mol) > 150: mL.append(Mol)
         mols.extend(mL)
     
     
     MyMolList.extend(mols)
     
     uniqps = {}
     for p in MyMolList:
        smi = Chem.MolToSmiles(Chem.AddHs(p))
        uniqps[smi] = p
     MyMolList = uniqps.values()
     
     print " The number of molecules is now: ", len(MyMolList)
     
     ####################################################
     uniqps = {}
     for p in MyMolList:
        smi = Chem.MolToSmiles(Chem.AddHs(p))
        uniqps[smi] = p
     MyMolList = uniqps.values()
     
     
     print " the number of compounds generated from the compound used to query: ", len(MyMolList)
     print Chem.MolToSmiles(MyMol)
     output_table=DataFrame()
    
     print len(MyMolList)

     for i,Molie in enumerate(MyMolList):
         print Chem.MolToSmiles(Molie)    
         output_table.set_value(i, 'SMILES',  Chem.MolToSmiles(Molie)) 


     from tabulate import tabulate
     #print tabulate(output_table, headers='keys', tablefmt='psql')


