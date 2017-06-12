#/bin/bash

for f in *.png; do echo "![alt tag](docs/$f?raw=true)" >> ../Readme.md ; done
