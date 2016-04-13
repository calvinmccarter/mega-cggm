#!/bin/sh
(cd ../Pseudo && make)
../Pseudo/pseudo_run -v 1 -y 0.1 -x 0.2 10 15 10 12 Yfile Xfile Lambdafile Thetafile statsfile
