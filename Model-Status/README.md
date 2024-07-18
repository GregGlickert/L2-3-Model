# L2-3-Model

This folder will contain all the updated write ups for the project for now. At some point when Nair gets to looking at them more we will have to convert into work docs. For now i think it is better and easier to keep in either MD or notebook format. We can then use pandoc to convert into docx.


## ways to convert notebook into better format
This will convert into html without showing python code better to do this than straight to docx
jupyter nbconvert connectivity-notebook.ipynb --no-input --to html
Then this will convert to word 
pandoc connectivity-notebook.html -s -o connectivity.docx