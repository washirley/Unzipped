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

#TestMol = Chem.MolFromSmiles(flow_variables['CompoundSmiles'])
#TestMol = Chem.MolFromSmiles('CC1OC(=NC1C1=NC(CS1)C1=NC(CS1)C(=O)O)c1ccccc1O')
# Kirromycin
#TestMol = Chem.MolFromSmiles('CCC(C(=O)NCC=CC=C(C)C(C(C)C1C(C(C(O1)C=CC=CC=C(C)C(=C2C(=O)C=CNC2=O)O)O)O)OC)C3(C(C(C(C(O3)C=CC=CC)(C)C)O)O)O')
# Barbamide: none current
#TestMol = Chem.MolFromSmiles('CC(C/C(=C\C(=O)N(C)C(CC1=CC=CC=C1)C2=NC=CS2)/OC)C(Cl)(Cl)Cl')
# Coelibactin
#TestMol = Chem.MolFromSmiles('CC1OC(=NC1C1=NC(CS1)C1=NC(CS1)C(=O)O)c1ccccc1O')
# Scabichelin
#TestMol = Chem.MolFromSmiles('NCCCC(C(=O)NC1CCCN(C1=O)O)NC(=O)C(CCCN(C(=O)C(NC(=O)C(CCCN(C(=O)C)O)NC)CO)O)NC')
# Benzene
#TestMol = Chem.MolFromSmiles('c1ccccc1')



#-----------------------
#FK-506
#TestMol = Chem.MolFromSmiles('CCC1C=C(CC(CC(C2C(CC(C(O2)(C(=O)C(=O)N3CCCCC3C(=O)OC(C(C(CC1=O)O)C)C(=CC4CCC(C(C4)OC)O)C)O)C)OC)OC)C)C')

#Barbamide
#TestMol = Chem.MolFromSmiles('CC(C/C(=C\C(=O)N(C)C(CC1=CC=CC=C1)C2=NC=CS2)/OC)C(Cl)(Cl)Cl')


#Kyrromycin
#TestMol = Chem.MolFromSmiles('CCC(C(=O)NCC=CC=C(C)C(C(C)C1C(C(C(O1)C=CC=CC=C(C)C(=C2C(=O)C=CNC2=O)O)O)O)OC)C3(C(C(C(C(O3)C=CC=CC)(C)C)O)O)O')

#-----------------------
TestMols=[]

#import Prism_Include
#exec(compile(open('Prism_Include.py').read()))
#execfile('Prism_Include.py')
execfile('Antibiotics_Include.py')
'''
TestMols.append(("FK-506",         Chem.MolFromSmiles('CCC1C=C(CC(CC(C2C(CC(C(O2)(C(=O)C(=O)N3CCCCC3C(=O)OC(C(C(CC1=O)O)C)C(=CC4CCC(C(C4)OC)O)C)O)C)OC)OC)C)C')))
TestMols.append(("Barbamide",      Chem.MolFromSmiles('CC(C/C(=C\C(=O)N(C)C(CC1=CC=CC=C1)C2=NC=CS2)/OC)C(Cl)(Cl)Cl')))
TestMols.append(("Kyrromycin",     Chem.MolFromSmiles('CCC(C(=O)NCC=CC=C(C)C(C(C)C1C(C(C(O1)C=CC=CC=C(C)C(=C2C(=O)C=CNC2=O)O)O)O)OC)C3(C(C(C(C(O3)C=CC=CC)(C)C)O)O)O')))
TestMols.append(("Tomamycin",      Chem.MolFromSmiles('CC=C1CC2C(NC3=CC(=C(C=C3C(=O)N2C1)OC)O)OC')))
TestMols.append(("Pyoluteorin",    Chem.MolFromSmiles('C1=CC(=C(C(=C1)O)C(=O)C2=CC(=C(N2)Cl)Cl)O')))
TestMols.append(("Fumitremorgin C",Chem.MolFromSmiles('CC(=CC1C2=C(CC3N1C(=O)C4CCCN4C3=O)C5=C(N2)C=C(C=C5)OC)C')))
TestMols.append(("Teicoplanin",    Chem.MolFromSmiles('CCCCCCCCCC(=O)NC1C(C(C(OC1OC2=C3C=C4C=C2OC5=C(C=C(C=C5)C(C6C(=O)NC(C7=CC(=CC(=C7C8=C(C=CC(=C8)C(C(=O)N6)NC(=O)C4NC(=O)C9C1=CC(=CC(=C1)O)OC1=C(C=CC(=C1)C(C(=O)NC(CC1=CC(=C(O3)C=C1)Cl)C(=O)N9)N)O)O)OC1C(C(C(C(O1)CO)O)O)O)O)C(=O)O)OC1C(C(C(C(O1)CO)O)O)NC(=O)C)Cl)CO)O)O')))
TestMols.append(("Ansamycin",      Chem.MolFromSmiles('CC1C=CC=C(C(=O)NC2=C3C(=NC4(N3)CCN(CC4)CC(C)C)C5=C6C(=C(C(=C5C2=O)O)C)OC(C6=O)(OC=CC(C(C(C(C(C(C1O)C)O)C)OC(=O)C)C)OC)C)C')))
TestMols.append(("Cerulenin",      Chem.MolFromSmiles('CC=CCC=CCCC(=O)C1C(O1)C(=O)N')))
TestMols.append(("Scabichelin",    Chem.MolFromSmiles('NCCCC(C(=O)NC1CCCN(C1=O)O)NC(=O)C(CCCN(C(=O)C(NC(=O)C(CCCN(C(=O)C)O)NC)CO)O)NC')))
TestMols.append(("Erythromycin",   Chem.MolFromSmiles('CCC1C(C(C(C(=O)C(CC(C(C(C(C(C(=O)O1)C)OC2CC(C(C(O2)C)O)(C)OC)C)OC3C(C(CC(O3)C)N(C)C)O)(C)O)C)C)O)(C)O')))
#
TestMols.append(("Tropane",  Chem.MolFromSmiles('CN1[C@@H]2CC[C@H]1CCC2')))
TestMols.append(("platensimide A",  Chem.MolFromSmiles('CC(=O)NCC[C@H](NC(=O)CC[C@]3(C)C(=O)C=CC24CC1CC(OC1(C)C2)C34)C(O)=O')))
TestMols.append(("Platensin",  Chem.MolFromSmiles('C[C@@]3(CCC(O)=O)C(=O)C=CC24CC1CC(OC1(C)C2)C34')))
TestMols.append(("Pironetin",  Chem.MolFromSmiles('C/C=C/C[C@H](C)[C@@H](OC)[C@@H](C)[C@H](O)C[C@H]1OC(=O)C=C(C)[C@H]1CC')))
TestMols.append(("Apicularin A",  Chem.MolFromSmiles('CC/C=C\C=C/C(=O)N\C=C\C[C@H]3C[C@H]1C[C@H](O)C[C@@H](O1)c2cccc(O)c2C(=O)O3')))
TestMols.append(("7-ACA",  Chem.MolFromSmiles('CC(=O)OCC2CS[C@@H]1[C@H](N)C(=O)N1C=2C(O)=O')))
TestMols.append(("Prostaglandin 2-alpha",  Chem.MolFromSmiles('CCCCC[C@H](O)/C=C/[C@H]1[C@H](O)C[C@H](O)[C@@H]1C\C=C/CCCC(O)=O')))
TestMols.append(("galanthamine",  Chem.MolFromSmiles('COc2ccc3CN(C)CC[C@@]14C=C[C@H](O)C[C@@H]1Oc2c34')))
TestMols.append(("Discodermolide",  Chem.MolFromSmiles('C[C@@H](C/C(/C)=C\[C@H](C)[C@@H](O)[C@@H](C)/C=C\[C@@H](O)C[C@@H]1OC(=O)[C@H](C)[C@@H](O)[C@H]1C)[C@@H](O)[C@H](C)[C@@H](OC(N)=O)[C@@H](C)/C=C\C=C')))
TestMols.append(("Halaven",  Chem.MolFromSmiles('CO[C@H]3[C@@H](C[C@H](O)CN)O[C@H]2C[C@H]9O[C@@H](CC[C@@H]8O[C@@H](CC[C@@]17C[C@H]6O[C@H]5[C@@H](O1)[C@H]4O[C@@H](CC(=O)C[C@@H]23)CC[C@@H]4O[C@H]5[C@H]6O7)CC8=C)C[C@@H](C)C9=C')))
TestMols.append(("shikimic acid",  Chem.MolFromSmiles('OC(=O)C1C[C@@H](O)[C@H](O)[C@H](O)C=1')))
TestMols.append(("salicylic acid",  Chem.MolFromSmiles('OC(=O)c1ccccc1O')))
TestMols.append(("morphine",  Chem.MolFromSmiles('CN3CC[C@@]15[C@H]2C=C[C@H](O)[C@@H]1Oc4c(O)ccc(C[C@H]23)c45')))
TestMols.append(("Quinine",  Chem.MolFromSmiles('C=C[C@H]1CN2CC[C@H]1C[C@H]2[C@H](O)c3ccnc4ccc(cc34)OC')))
#TestMols.append(("Indigo",  Chem.MolFromSmiles('O=C1/C(/Nc2ccccc12)=C4\Nc3ccccc3C\4=O')))
TestMols.append(("Indigo",  Chem.MolFromSmiles('O=C1\C(NC2=C1C=CC=C2)=C3/NC4=C(C=CC=C4)C3=O')))
TestMols.append(("Neomycin B",  Chem.MolFromSmiles('NC[C@@H]4O[C@H](O[C@H]3[C@@H](O)[C@H](O[C@@H]1[C@@H](O)[C@H](N)C[C@H](N)[C@H]1O[C@H]2O[C@H](CN)[C@@H](O)[C@H](O)[C@H]2N)O[C@@H]3CO)[C@H](N)[C@@H](O)[C@@H]4O')))
TestMols.append(("artemisinin",  Chem.MolFromSmiles('C[C@@H]1CC[C@H]2[C@@H](C)C(=O)O[C@@H]3O[C@@]4(C)CC[C@@H]1[C@@]23OO4')))
TestMols.append(("cholesterol",  Chem.MolFromSmiles('CC(C)CCC[C@@H](C)[C@H]3CC[C@H]4[C@@H]2CC=C1C[C@H](O)CC[C@]1(C)[C@H]2CC[C@]34C')))
TestMols.append(("coumarin",  Chem.MolFromSmiles('O=C2C=Cc1ccccc1O2')))
TestMols.append(("nicotine",  Chem.MolFromSmiles('CN1CCC[C@H]1c2cccnc2')))
TestMols.append(("platencin subst",  Chem.MolFromSmiles('C[C@@]1(CCC(N)=O)C(=O)C=C[C@]23CC[C@@H](C[C@@H]12)C(=C)C3')))
TestMols.append(("Vincadifformine core",  Chem.MolFromSmiles('COC(=O)C2C[C@@H]4CCN5CC[C@@]3(c1ccccc1NC=23)[C@H]45')))
TestMols.append(("Penicillin core",  Chem.MolFromSmiles('CC(=O)N[C@@H]1C(=O)N2[C@@H](C(O)=O)C(C)(C)S[C@H]12')))
'''
#import Prism_Include


for TestMolName, TestMol in TestMols:
     print "Starting Mol: ", TestMolName
     # determine the cut-offs based on the MW and Number of Heavy Atoms
     MolWt = Descriptors.MolWt(TestMol)
     HeavyAtoms = Descriptors.HeavyAtomCount(TestMol) 
     
     weights_Param_1 = [ 0.00198015, -0.02093066]
     weights_Param_2 = [ 0.01696928, -0.12057078]
     Desc = [MolWt, HeavyAtoms]
     Param_1 = np.dot(Desc, weights_Param_1) + 0.35198505514161238
     Param_2 = np.dot(Desc, weights_Param_2) + 4.1740354586675785
     
     MyMol=Chem.RemoveHs(TestMol)
     R=MyMol.GetRingInfo()
     NumRings = R.NumRings()
     
     # useHs=False,branchedPaths=False
     MyMolFp = Chem.RDKFingerprint(MyMol,fpSize=2048, minPath=7, useBondOrder=False, nBitsPerHash=1, minSize=2048, useHs=False)
     #MyMolFp = Chem.PatternFingerprint(MyMol,fpSize=2048)
     MyMolList = [MyMol]
     #####################################################################################################
     ##########################  Create the AntiSmash List  ###############################################
     
     supp = []
     
     List = ['test.smi','new_public_microbes.smi','new_npugenomes.smi']
     
     suppl1 = Chem.SmilesMolSupplier(List[1], smilesColumn=1, nameColumn=0,titleLine=True,sanitize=False)
     suppl2 = Chem.SmilesMolSupplier(List[2], smilesColumn=1, nameColumn=0,titleLine=True,sanitize=False)
     
     for me in suppl1:
         me.GetNumAtoms()
         Chem.SanitizeMol(me, (Chem.SanitizeFlags.SANITIZE_ADJUSTHS | Chem.SanitizeFlags.SANITIZE_SYMMRINGS))
         me = Chem.RemoveHs(me)
         if Descriptors.MolWt(me) > 125: supp.append(me)
     
     for me in suppl2:
         me.GetNumAtoms()
         Chem.SanitizeMol(me, (Chem.SanitizeFlags.SANITIZE_ADJUSTHS | Chem.SanitizeFlags.SANITIZE_SYMMRINGS))
         me = Chem.RemoveHs(me)
         if Descriptors.MolWt(me) > 125: supp.append(me)
         
     print " Total number of Fltered AntiSmash predictions: ", len(supp)
     
     ######################## Fingerprint the AntiSMASH structures ##########################
     # useHs=False, branchedPaths=False
     fps = [Chem.RDKFingerprint(x,fpSize=2048, minPath=7, useBondOrder=False, nBitsPerHash=1, minSize=2048, useHs=False) for x in supp]
     
     #fps = [Chem.PatternFingerprint(x,fpSize=2048) for x in supp]
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
     #        Chem.Kekulize(AMol)
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
     #                m = Chem.RemoveHs(m)
                     AList.append(m)
         MyMolList.extend(AList)
     
     uniqps = {}
     for p in MyMolList:
     #   smi = Chem.MolToSmiles(Chem.RemoveHs(Chem.AddHs(p)))
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
     #       AllChem.SetAromaticity(m)
     #       m = Chem.AddHs(m)
     #       m = Chem.RemoveHs(m)
                Chem.SanitizeMol(m, (Chem.SanitizeFlags.SANITIZE_ADJUSTHS | Chem.SanitizeFlags.SANITIZE_SYMMRINGS))
                if Descriptors.MolWt(m) > 150: mL.append(m)
             mols.extend(mL)
     
         m = AllChem.DeleteSubstructs(AMol,pattD)
         Chem.MolToSmiles(m,kekuleSmiles=True)
     #    Chem.Kekulize(m)
         AllChem.Cleanup(m)
         AllChem.SetAromaticity(m)
     #    m = Chem.RemoveHs(m)
         mols.append(m)
     
         mL = []
         m = AllChem.ReplaceSubstructs(AMol,pattR,replR)
         for Rot,Mol in enumerate(m):
             Mol.UpdatePropertyCache(strict=False)
             AllChem.Cleanup(Mol)
             AllChem.SetAromaticity(Mol)
     #        Mol = Chem.RemoveHs(Mol)  
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
     # useHs=False,branchedPaths=False
     fpq = [Chem.RDKFingerprint(x,fpSize=2048, minPath=7,  useBondOrder=False,nBitsPerHash=1, minSize=2048,useHs=False) for x in MyMolList]
     
     #fpq = [Chem.PatternFingerprint(x,fpSize=2048) for x in MyMolList]
     AllQBits = [list(fpqe.GetOnBits()) for fpqe in fpq]
     AllQBits = list(set([j for i in AllQBits for j in i]))
     if len(AllQBits) == 0:
         AllQBits=[2047]
     bv = DataStructs.ExplicitBitVect(2048)
     bv.SetBitsFromList(AllQBits)
         
     
     print " the number of compounds generated from the compound used to query: ", len(MyMolList)
     print len(AllQBits)
     print Chem.MolToSmiles(MyMol)
     print len(MyMolFp.GetOnBits())
     
     ################################  End Query Ensemble  Creation ################################################################
     ############################  Find the AntiSmash compounds with the most overlap with the query ensemble fingerprint  #########
     
     print " Culling the AntiSmash compounds "
     print " examining fingerprints : ", len(fps)
     MySearchList=[]
     similarity = 0.17
     print similarity, (len(MySearchList) == 0)
     # am testing the while list to pick up compounds specficially for s.; also am adjusitng the < 1.5 times the bit lengths; this was 1
     while len(MySearchList) == 0:
         print " in the list generation loop "
         
         for i,fp in enumerate(fps):
             andfp = bv&fp
             obl = list(andfp.GetOnBits())
     #    if len(obl) > 0.90*len(AllQBits):
     #    if len(obl) > 0.70*len(AllQBits):
     #    if len(obl) > 0.80*len(AllQBits):
     #    if len(obl) > 0.85*len(AllQBits):
     ##        print " test x ", len(obl), 0.25*len(AllQBits), len(fp.GetOnBits()), len(AllQBits),  DataStructs.FingerprintSimilarity(bv, fp)
             if len(obl) > 0.25*len(AllQBits):
                if len(fp.GetOnBits()) < 1.5*len(AllQBits):
     #           print DataStructs.FingerprintSimilarity(bv, fp),  DataStructs.FingerprintSimilarity(bv, fp, metric=DataStructs.DiceSimilarity)
                    if  DataStructs.FingerprintSimilarity(bv, fp) > similarity:
                         MySearchList.append(supp[i])
     #                    print " test x ", len(obl), 0.25*len(AllQBits), len(fp.GetOnBits()), len(AllQBits),  DataStructs.FingerprintSimilarity(bv, fp)
     #    MySearchList.append(supp[i])
         similarity = similarity - 0.005
         print similarity
         
         if similarity <= 0.01:
             MySearchList.append(supp[0])
             break
     
     
     
     print " length: ", len(MySearchList), Chem.MolToSmiles(MySearchList[0])
     for x in MySearchList:  print Chem.MolToSmiles(x),"\n" 
     
     # ,branchedPaths=False, useHs=FalseFalse
     fps = [Chem.RDKFingerprint(x,fpSize=2048, minPath=7, useBondOrder=False,nBitsPerHash=1, minSize=2048, useHs=False) for x in MySearchList]
     #fps = [Chem.PatternFingerprint(x,fpSize=2048) for x in MySearchList]
     
     print len((MyMolFp&fp).GetOnBits()), len(obl)
     print " Number of AntiSmash compounds to further search: ", len(MySearchList)
     ######################################## Reduce the Ensemble of Query Compounds  ###############################################
     
     MySkinnyMolList = []
     scale = 0.25
     while len(MySkinnyMolList) == 0:
         for number, fp in enumerate(fpq):
           for fs in fps:
     #    if DataStructs.FingerprintSimilarity(fp,fs, metric=DataStructs.TanimotoSimilarity) > 0.90*Param_1: 
             if DataStructs.FingerprintSimilarity(fp,fs, metric=DataStructs.TanimotoSimilarity) > scale*Param_1: 
                MySkinnyMolList.append(MyMolList[number])
                break
         scale += -0.05
         
     print MyMolList[1]
     print len(MySkinnyMolList)
     print len(MyMolList)
     
     
     uniqpss = {}
     for p in MySkinnyMolList:
        smi = Chem.MolToSmiles(p)
        uniqpss[smi] = p
     MyMolList = uniqpss.values()
     
     print " now about to start the MCS with: ", len(uniqpss), " of the ensemble, against AntiSmash: ", len(MySearchList)
     
     ########################################  Below is the MCS section  #######################################################
     i=0
     output_table=DataFrame({'matches': ["Data"]}, columns= ['matches'])
     for j,p in enumerate(MySearchList):
       p.UpdatePropertyCache(strict=False)
       for q in MyMolList:
           q.UpdatePropertyCache(strict=False)
     
           MolTestSet=[p,q]
     #      res=rdFMCS.FindMCS(MolTestSet, bondCompare=rdFMCS.BondCompare.CompareAny, atomCompare=rdFMCS.AtomCompare.CompareAny,timeout=10)
           res=rdFMCS.FindMCS(MolTestSet, bondCompare=rdFMCS.BondCompare.CompareAny, atomCompare=rdFMCS.AtomCompare.CompareElements,timeout=1)
           if res.numAtoms >= 1.2*Param_2:  
     
               query = Chem.MolFromSmarts(res.smartsString)
               highlightAtoms=list(p.GetSubstructMatch(query))
               output_table.set_value(i, 'TestMolName', TestMolName)
               output_table.set_value(i, 'matches',  p.GetProp('_Name'))
               output_table.set_value(i, 'SMILES',  Chem.MolToSmiles(p))
               output_table.set_value(i, 'MCS Atoms', str(highlightAtoms))
               output_table.set_value(i, 'MCS Smarts', res.smartsString)
               i += 1
               break
     ########################################  End of  the MCS section  #######################################################
     
     from tabulate import tabulate
     print tabulate(output_table, headers='keys', tablefmt='psql')
     print "Ending Mol: ", TestMolName
