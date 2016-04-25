todo list 

things that could be done (in random order):
1. add function to automate selection of helix points.[footnote1] from data.[footnote2] 
2. 


footnotes (on order of appearance):
#-----------------------------------------------------------------------------------------------------------------------------------------
*1. okay, what do i mean by automate? great question. i have some ideas about how to do it, but the sake of real progress, i-ll not share
here (yet?) because they-re not as awesome as i-d like (super empirical at this point). that very well may be the only solution-type, 
however. i mean, for different problems, such as deriving helical parameters of single stranded RNA (ssRNA) in the ribosome (goto rcsb.org
and search "5HCQ" -- if you-re really bout-it bout-it, install VMD from ks.uiuc.edu/* and visualize the ribosome interactively), we may
need to develop a set of RNA-specific selection rules. actually, 99% sure that-s exactly what needs doing based on the literature...


#-----------------------------------------------------------------------------------------------------------------------------------------
*2. there are multiple formats the data can be read, or even streamed, that needs code. 
formats that helios currently "supports"

(1) "PDB" Protein Data Bank (see http://deposit.rcsb.org/adit/docs/pdb_atom_format.html --- namely, the "ATOM" section. fortran ftw? sigh.)
code: Kreadpdb.f90
input parameter to main code that toggles this format: "coord_type = 1"

(2) "xyz" rudest format in which order triples (Python calls this a list comprehension, i think --- see http://stackoverflow.com/questions/1403674/pythonic-way-to-return-list-of-every-nth-item-in-a-larger-list)
the python version of the code would look something like:
with open ('filename') as input:
  for line in input:
    coordinates = [ item.strip() for item in line.split() ]
#and then i think the x y z --- assuming that-s how they-re ordered in 'filename'! --- gets broken up like this *i think*:
x = coordinates[0:3] #starting with first element in array (or "list"), get every third thenceforth
y = coordinates[1:3] #starting with second element in array, get every third thence..
z = coordinates[2:3] #starrting with third element in array, get every third thence...

#-----------------------------------------------------------------------------------------------------------------------------------------
*3.....
