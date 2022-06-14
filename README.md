**The only file by me is 'squeeze.py', the rest comes from https://gitea.tforgione.fr/tforgione/obja**

The program is an implementation of the SQUEEZE mesh compression method described here : https://www.researchgate.net/publication/3852440_SQUEEZE_fast_and_progressive_decompression_of_triangle_meshes  


## Creation of obja file
- Launch squeeze.py with the obj file name as argument. ex : "py squeeze.py suzanne" for the models/suzanne.obj file


## Visualisation
- Launch the server : py obja/server.py
- Open web page : localhost:8000/?output/suzanne.obja


## Change display speed
- Open js/main.js and change the last value line 10 : 'loader = new Loader(url, 1024, 20);' 
