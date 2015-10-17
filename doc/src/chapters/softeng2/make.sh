#!/bin/sh
name=softeng2

# Compile locally
#bash ../make.sh $name

# Compile multiple formats and publish documents and soure code to INF5620
bash ../make.sh $name sphinx publish src
