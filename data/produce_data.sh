#/bin/bash

# Runs the algorithm on all of the matrices and saves the output.

declare -rx DEVROOT=~/scratch/sgd
declare -rx LIBROOT="${DEVROOT}/lib"
declare -rx REPOROOT="${DEVROOT}/SGD"
declare -rx SRCROOT="${REPOROOT}/cpp"
declare -rx MATRICES_ROOT="${REPOROOT}/matrices"
declare -rx DATAROOT="${REPOROOT}/data"

declare -rx MATRIX_LIST="${MATRICES_ROOT}/matrix_list.txt"
declare -rxa THREAD_COUNTS=(1 2 4 8 16)

declare -rx PRNG_SEED=0

# How many major iterations (sequences of n steps) to do.
declare -rx MI_CNT=1000

# How many major iterations are there per epoch. The residual is computed after
# each epoch.
declare -rx MIS_PER_EPOCH=1

# Pattern of filenames: <directory>/<matrix>-T<threads>.txt
declare -rx PATTERN="%s/%s-T%02d.txt"

driver() {
    build
    cat "${MATRIX_LIST}" |run
}

# Builds the code.
build() {
    pushd "${PWD}"
    cd "${SRCROOT}"
    CPATH="${LIBROOT}/pfunc/include":"${LIBROOT}/boost" make
    popd
}

# Runs the code on the matrices it reads from stdin.
run() {
    while read M; do
        for T in "${THREAD_COUNTS[@]}"; do
            local -x MATRIX="${MATRICES_ROOT}/${M##*/}.mtx"
            local -x RHS="${MATRICES_ROOT}/y-${M##*/}.txt"
            local -x DATAFILE=`printf "${PATTERN}" "${DATAROOT}" "${M##*/}" "${T}"`
            local -x SIZE=`cat "${MATRIX}" | matrix_size`
            local -x EPOCHS=$(( MI_CNT / MIS_PER_EPOCH ))
            local -x EPOCH_SIZE=$(( MIS_PER_EPOCH * SIZE ))
            "${SRCROOT}/harness-opt" \
                seed "${PRNG_SEED}" \
                num-threads "${T}" \
                max-epochs "${EPOCHS}" \
                epoch-size "${EPOCH_SIZE}" \
                A-file-path "${MATRIX}" \
                Y-file-path "${RHS}" \
                M "${SIZE}" \
                | tail --lines=+2 > "${DATAFILE}"
        done
    done
}

# Prints the number of rows in a matrix. Reads the matrix in Matrix Market
# program format from stdin.
matrix_size() {
    awk '/^[^%]/ { print $1; exit }'
}

driver
