#! /bin/bash

echo SPAWN "$@"

export MV2_ENABLE_AFFINITY=0
export DYLD_LIBRARY_PATH="$VISTLE_DYLD_LIBRARY_PATH"

if [ -n "$SLURM_JOB_ID" ]; then
   exec mpiexec -bind-to none "$@"
   #exec srun --overcommit "$@"
   #exec srun --overcommit --cpu_bind=no "$@"
fi

case $(hostname) in
   viscluster70)
      if [ -z "$MPIHOSTS" ]; then
         MPIHOSTS=$(echo viscluster70 viscluster{51..60} viscluster{71..77}|sed -e 's/ /,/g')
      fi
      if [ "$MPISIZE" = "" ]; then
         MPISIZE=2
      fi
      #exec mpirun -np ${MPISIZE} -hosts ${MPIHOSTS} "$@"
      exec mpirun -np ${MPISIZE} -hosts ${MPIHOSTS} -bind-to none -envall "$@"
      #exec mpirun -np 8  -hosts localhost "$@"
      ;;
   viscluster*)
      if [ -z "$MPIHOSTS" ]; then
         MPIHOSTS=$(echo viscluster{50..60} viscluster{71..75}|sed -e 's/ /,/g')
         #MPIHOSTS=$(echo viscluster{50..60}|sed -e 's/ /,/g')
      fi
      if [ "$MPISIZE" = "" ]; then
         MPISIZE=16
      fi
      #echo LD_LIBRARY_PATH=$LD_LIBRARY_PATH
      #exec mpirun -np ${MPISIZE} -hosts ${MPIHOSTS} -bind-to none -envlist LD_LIBRARY_PATH "$@"
      exec mpirun -np ${MPISIZE} -hosts ${MPIHOSTS} -bind-to none -envall "$@"
      ;;
   *)
      #echo mpirun "$@"
      #exec mpirun -np 1 "$@"
      #echo EXEC: "$@"
      #exec mpirun -np 2 "$@"
      #exec xterm -e gdb --args "$@"
      if [ "$MPISIZE" = "" ]; then
         MPISIZE=1
      fi
      if [ -z "$MPIHOSTS" ]; then
         echo mpirun -np ${MPISIZE} "$@"
         exec mpirun -envall -np ${MPISIZE} "$@"
      else
         echo mpirun -np ${MPISIZE} -hosts ${MPIHOSTS} "$@"
         exec mpirun -envall -np ${MPISIZE} -hosts ${MPIHOSTS} "$@"
      fi
      ;;
esac

exec "$@"
