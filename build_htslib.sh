tar -xjf htslib-1.16.tar.bz2
cd htslib-1.16
autoheader
autoconf
./configure --prefix=`pwd`
make
make install
