#!/bin/bash

# Downloads matrices from UFL and creates RHS files for them.

declare -rx MATLAB="/gsa/yktgsa/projects/w/watapps/matlab/linux/R2009a/bin/matlab"
declare -rx UFGET="http://www.cise.ufl.edu/research/sparse/mat/UFget.tar.gz"
declare -rx MATRIX_LIST=matrix_list.txt
declare -rx URL_PATTERN="http://www.cise.ufl.edu/research/sparse/MM/%s.tar.gz"

driver() {
    cat "${MATRIX_LIST}" | download_mm
    cat "${MATRIX_LIST}" | make_rhs
}

# Downloads matrices in Matrix Market format. Reads matrix names from stdin.
download_mm() {
    while read M; do
        url=`printf "${URL_PATTERN}" "$M"`
        file="${M##*/}".tar.gz
        wget "${url}" && \
            tar x --strip-components=1 --file "${file}" && \
            rm "${file}"
    done
}

# Uses Matlab to produce right-hand side vectors.
make_rhs() {
    wget "${UFGET}" && tar xf "${UFGET##*/}" && rm "${UFGET##*/}"
    fmt_matrix_list | \
        make_rhs_script | \
        "${MATLAB}" -nodisplay
    rm -rf UFget
}

# Produces a Matlab script that defines a cell array of matrix names. Reads matrix
# names from stdin.
fmt_matrix_list() {
    printf "MATRICES = { "
    if read M; then
        printf "'${M}'"
    fi
    while read M; do
        printf ", '${M}'"
    done
    printf " };\n"
}

# Produces a Matlab script for generating right-hand side vectors. Reads a cell
# array of matrix names from stdin.
make_rhs_script() {
    cat <<EOF
addpath('UFget');
EOF

    cat -

    cat <<EOF
PRECISION = '%30.16e';
for M = MATRICES
  A = UFget(M{1});
  [ig, name, ig, ig] = fileparts(M{1});
  rhs = ['y-' name '.txt'];
  dlmwrite(rhs, A.A * ones(size(A.A, 1), 1), 'precision', PRECISION);
end
EOF
}

driver
