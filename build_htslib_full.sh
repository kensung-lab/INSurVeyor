tar -xjf htslib-1.19.1.tar.bz2
cd htslib-1.19.1
./configure --prefix=`pwd` --enable-libcurl
make
make install
