tar -xjf htslib-1.13.tar.bz2
cd htslib-1.13
autoheader
autoconf
./configure --prefix=`pwd` --disable-bz2 --disable-lzma
make
make install
