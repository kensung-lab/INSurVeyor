tar -xjf htslib-1.20.tar.bz2
cd htslib-1.20
./configure --prefix=`pwd` --enable-libcurl
make
make install
