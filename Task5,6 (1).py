#!/usr/bin/env python
# coding: utf-8

# In[1]:


from pyopenms import * 


# In[2]:


# Proteolytic Digestion with Trypsin yeast
dig = ProteaseDigestion()
dig.getEnzymeName()
bsa = "".join([l.strip() for l in open("uniprot-yourlist_M202111205BF3C56A578D7D6DFD1FC81EE5DA77302FCEBC7.fasta").readlines()[1:]])
bsa = AASequence.fromString(bsa)
# create all digestion products
result = []
dig.digest(bsa, result)
for s in result:
    print(s.toString())
len(result) 


# In[3]:


# Proteolytic Digestion with Trypsin Human 
dig = ProteaseDigestion()
dig.getEnzymeName()
bsa = "".join([l.strip() for l in open("uniprot-yourlist_M202112026320BA52A5CE8FCD097CB85A53697A35338CFBY.fasta").readlines()[1:]])
bsa = AASequence.fromString(bsa)
# create all digestion products
result = []
dig.digest(bsa, result)
for s in result:
    print(s.toString())
len(result) 


# In[5]:


# Y-ion spectrum 
tsg = TheoreticalSpectrumGenerator()
spec1 = MSSpectrum()
peptide = AASequence.fromString("PEPTIDERDLQMTQSPSSLSVSVGDRPEPTIDE")
# standard behavior is adding b- and y-ions of charge 1
p = Param()
p.setValue("add_b_ions", "false")
p.setValue("add_metainfo", "true")
tsg.setParameters(p)
tsg.getSpectrum(spec1, peptide, 1, 1) # charge range 1:1

# Iterate over annotated ions and their masses
print("Spectrum 1 of", peptide, "has", spec1.size(), "peaks.")
for ion, peak in zip(spec1.getStringDataArrays()[0], spec1):
    print(ion.decode(), "is generated at m/z", peak.getMZ())


# In[ ]:


# Y-ion spectrum yeast
dig = ProteaseDigestion()
dig.getEnzymeName()
bsa = "".join([l.strip() for l in open("uniprot-yourlist_M202111205BF3C56A578D7D6DFD1FC81EE5DA77302FCEBC7.fasta").readlines()[1:]])
bsa = AASequence.fromString(bsa)
# create all digestion products
result = []
dig.digest(bsa, result)
peptides = [AASequence.fromString(s.toString()) for s in result]

# Iterate over annotated ions and their masses
for peptide in peptides:
    tsg = TheoreticalSpectrumGenerator()
    spec1 = MSSpectrum()

    # standard behavior is adding b- and y-ions of charge 1
    p = Param()
    p.setValue("add_b_ions", "false")
    p.setValue("add_metainfo", "true")
    tsg.setParameters(p)
    tsg.getSpectrum(spec1, peptide, 1, 1) # charge range 1:1
    print("Spectrum 1 of", peptide, "has", spec1.size(), "peaks.")
    for ion, peak in zip(spec1.getStringDataArrays()[0], spec1):
        print(ion.decode(), "is generated at m/z", peak.getMZ())


# In[ ]:




