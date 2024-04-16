#!/usr/bin/env bash

version=4.2.7

docker build -t bcl_convert-$version .
docker tag bcl_convert-$version us-docker.pkg.dev/microbiome-xavier/broad-microbiome-xavier/bcl_convert:$version
docker push us-docker.pkg.dev/microbiome-xavier/broad-microbiome-xavier/bcl_convert:$version