docker run --rm -v ${pwd}:/data -p 8888:8888 --entrypoint=""  kavrakilab/hla-arena:sars-arena jupyter notebook --port=8888 --no-browser --ip=0.0.0.0 --allow-root
