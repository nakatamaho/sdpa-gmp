# sdpa-gmp
# how to build
I verified build on Ubuntu 20.04
```
rm -rf sdpa-gmp
git clone https://github.com/nakatamaho/sdpa-gmp.git
cd sdpa-gmp
aclocal ; autoconf ; automake --add-missing
autoreconf --force --install
./configure --enable-openmp=yes --enable-shared=yes
make -j4
```
