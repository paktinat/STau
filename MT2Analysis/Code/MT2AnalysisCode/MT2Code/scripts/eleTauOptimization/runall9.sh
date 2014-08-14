for s in 0 1 2 3 4
do
  ./STauEleTauCutOptimization.py -s $s -d /nfs-6/userdata/hbakhshi/  >& out$s_2  &
done
