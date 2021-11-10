#! /bin/sh
#
# entrypoint.sh
# Copyright (C) 2021 dlilien <dlilien90@gmail.com>
#
# Distributed under terms of the MIT license.
#


export PLAT=manylinux2014_x86_64

# Do a normal build
for PYBIN in /opt/python/cp3[7-9]-cp*/bin; do
    "${PYBIN}/pip" install numpy==1.19.0 cython;
    "${PYBIN}/pip" wheel --no-deps -w /github/workspace/wheelhouse/ .;
done

# Make the wheels into manylinux
ls wheelhouse/*.whl
for whl in wheelhouse/*.whl; do
    auditwheel repair "$whl" --plat $PLAT -w /github/workspace/dist/;
done
