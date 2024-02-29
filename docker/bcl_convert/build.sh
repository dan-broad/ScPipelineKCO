#!/usr/bin/env bash

version=4.2.7

docker build -t bcl_convert-$version .
docker tag bcl_convert-$version gcr.io/microbiome-xavier/bcl_convert:$version
docker push gcr.io/microbiome-xavier/bcl_convert:$version