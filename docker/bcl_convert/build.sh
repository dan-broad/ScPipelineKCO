#!/usr/bin/env bash

docker build -t bcl_convert-3.10.5 .
docker tag bcl_convert-3.10.5 gcr.io/microbiome-xavier/bcl_convert:3.10.5
docker push gcr.io/microbiome-xavier/bcl_convert:3.10.5
