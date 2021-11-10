FROM quay.io/pypa/manylinux2014_x86_64
COPY entrypoint.sh /entrypoint.sh
RUN /opt/python/cp37-cp37m/bin/pip install twine
ENTRYPOINT ["/entrypoint.sh"]
