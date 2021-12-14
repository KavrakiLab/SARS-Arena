# SARS-Arena
SARS-Arena: A Pipeline for Selection and Structural HLA Modeling of Conserved Peptides of SARS-relate

## Installation
1. If you don’t already have it, install Docker.

  Docker for Mac or Windows: https://www.docker.com/products/docker-desktop

  Docker for Linux: https://docs.docker.com/install

2. In a command prompt, pull the SARS-Arena image from Docker Hub by typing:

  docker pull kavrakilab/hla-arena:sars-arena

3. Run the following command in the terminal:

  `./start_SARS-Arena.sh`

4. This should generate a URL with the following format:

  `http://127.0.0.1:8888/?token=<token_value>`

4. Copy and paste this URL into a browser, and open any available Jupyter notebook (i.e., one of the files with extension .ipynb). Note that all the data created in the container will be saved inside sub-directories of the current folder.

5. Check out the [DOCUMENTATION](), also provided alongside the Jupyter notebooks, for additional information on the workflows and available functions.

Enjoy SARS-Arena!

 
