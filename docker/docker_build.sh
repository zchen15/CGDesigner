#!/usr/bin/env bash
git add .;git commit -m 'rebuild';git push
docker login -u="\$app" -p="TTFIWPVU2ZHF4VWAWEP4HHOP4HUQ9I9KT6P4EKASTO8TU1JWT6W88XTDJB9R56I0YA926QOULEBI5FA7IGBC1IFYKLGO12VITRB1CIBH42HXAZ3GBDPLME7P" quay.bopcluster.net
docker build -t quay.bopcluster.net/zchen15/embryobs .
docker push quay.bopcluster.net/zchen15/embryobs
