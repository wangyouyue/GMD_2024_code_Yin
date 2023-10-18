#!/bin/sh
JID=`pjsub -z jid --step --sparam "sn=1" -N $(basename $PWD) SDM_job.sh`
cd ./results
pjsub --step --sparam "jid=${JID%??}, sn=2, sd=ec!=0:after:1" ncl.sh
pjsub --step --sparam "jid=${JID%??}, sn=3, sd=ec!=0:after:2" merge.sh