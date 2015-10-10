#!/bin/sh
(cd ../AltNewtonBCD && make)
../AltNewtonBCD/hugecggm_run -y 0.1 -x 0.2 10 15 10 12 Yfile Xfile Lambdafile Thetafile statsfile
