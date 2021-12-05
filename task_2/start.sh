for nproc in 2 4 8 16 64
do
    for eps in 1e-4 2e-5 8e-6
    do
        mpisubmit.bg -n $nproc main --stdout res_$eps/$nproc.txt -- $eps
    done
done

