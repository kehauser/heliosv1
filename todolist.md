todo list 

things that could be done (in random order):
1. add function to automate selection of helix points.[footnote1] from data.[footnote2] 
2. add function to better handle file formats, including netcdf
3. add function to determine if a helix is kinked
4. add function to calculate the angle between the helical axes subtended by the kink between them
5. add GUI, i guess
6. add function to calculate bending.[footnote3]


footnotes (in order of appearance):
#--------------------------------------------------------------------------------------------------------------------------------
*1. okay, what do i mean by automate? great question. i have some ideas about how to do it, but the sake of real progress, i-ll not share
here (yet?) because they-re not as awesome as i-d like (super empirical at this point). that very well may be the only solution-type, 
however. i mean, for different problems, such as deriving helical parameters of single stranded RNA (ssRNA) in the ribosome (goto rcsb.org
and search "5HCQ" -- if you-re really bout-it bout-it, install VMD from ks.uiuc.edu/* and visualize the ribosome interactively), we may
need to develop a set of RNA-specific selection rules. actually, 99% sure that-s exactly what needs doing based on the literature...


#--------------------------------------------------------------------------------------------------------------------------------
*2. there are multiple formats the data can be read, or even streamed, that needs code. 
formats that helios currently "supports"

(1) "PDB" Protein Data Bank (see http://deposit.rcsb.org/adit/docs/pdb_atom_format.html --- namely, the "ATOM" section. 
fortran ftw? sigh.)
code: Kreadpdb.f90
input parameter to main code that toggles this format: "coord_type = 1"

(2) "xyz" rudest format in which order triples (Python calls this a list comprehension, i think --- see 
http://stackoverflow.com/questions/1403674/pythonic-way-to-return-list-of-every-nth-item-in-a-larger-list)
the python version of the code would look something like:
with open ('filename') as input:
  for line in input:
    coordinates = [ item.strip() for item in line.split() ]
#and then i think the x y z --- assuming that-s how they-re ordered in 'filename'! --- gets broken up like this:
x = coordinates[0:3] #starting with first element in array (or "list"), get every third thenceforth
y = coordinates[1:3] #starting with second element in array, get every third thence..
z = coordinates[2:3] #starting with third element in array, get every third thence...

#--------------------------------------------------------------------------------------------------------------------------------
*3. ok, so this totally depends on knowing exactly the behaviour of our method under varying irregularities in twist, rise, 
radius and number points/turn. depending on how that goes, let-s assume for protein helical secondary structure 
(wikipedia.org if this is new for you) we need at least 8 consecutive CA atoms to define the helix. in increments of one CA atom along 
the protein chain (assume this is >> 8 CA atoms) at a time, derive the helical axis and store its direction in an array. continue
this all the way down the helix (saving the angles for only those 8-CA atom helices that pass muster, having reasonable fit-residuals
and displaying helical parameters consistent with alpha-helix geometry). go back into the array and calculate the angle between the 
helical axis vectors of successive helical elements.

Questions? Comments? Insults?
email me at kevin.e.hauser.ATCG.gmail.coom






