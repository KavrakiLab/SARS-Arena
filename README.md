# SARS-Arena
SARS-Arena: Sequence and Structure-Guided Selection of Conserved Peptides from SARS-related Coronaviruses for Novel Vaccine Development

## Installation
1. If you don’t already have it, install Docker.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Docker for Mac or Windows: https://www.docker.com/products/docker-desktop

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Docker for Linux: https://docs.docker.com/install

&nbsp;

2. In a command prompt, pull the SARS-Arena image from Docker Hub by typing:

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;*docker pull kavrakilab/hla-arena:sars-arena*

&nbsp;

3. Clone this repo and run the following command in the terminal:

  `./start_SARS-Arena.sh`
  
  &nbsp;

4. This should generate a URL with the following format:

  `http://127.0.0.1:8888/?token=<token_value>`
  
  &nbsp;

5. Copy and paste this URL into a browser, and open the Jupyter notebook (i.e., one of the files with extension .ipynb). Note that all the data created in the container will be saved inside sub-directories of the current folder.

&nbsp;

6. Check out the [DOCUMENTATION](https://github.com/KavrakiLab/SARS-Arena/edit/master/DOCUMENTATION.html), also provided alongside the Jupyter notebooks, for additional information on the workflows and available functions.

&nbsp;

Enjoy SARS-Arena!

 
