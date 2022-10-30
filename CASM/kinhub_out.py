
from CASM.load_tree import get_kinase_name_map

arr = ['Q9UQM7', 'Q16566', 'P17612', 'Q13554', 'Q13153', 'P51812']	
arr2=['O43683', 'Q9H5K3', 'Q8IWB6', 'O60566', 'Q86YV5', 'Q96S44', 'Q8IV63', 'O43187', 'Q99986', 'Q8IZE3', 'Q86Y07', 'Q8TF76', 'Q7Z7A4']
	

h1 = """
# @ triangle : size=30 : fill=rgb(240,230,140) : stroke=red : strokeWidth = 3.0
@ 3 : 30 : 240,230,140 : red : 3
"""

h2 = """
# square : size=25 : fill=magenta : stroke=black (default) : strokeWidth = 1.5 (default)
@ 4 : 25 : magenta
"""


d = get_kinase_name_map(
    "datasets/kinase_table.txt",
    "UniprotID",
)

print(h1)
for i in arr:
    print(d[i]["xName"])

print()
print(h2)
for i in arr2:
    print(d[i]["xName"])