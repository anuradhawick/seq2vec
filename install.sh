echo "Downloading boost"
cwd=$(pwd)
wget --continue --timestamping --no-clobber https://boostorg.jfrog.io/artifactory/main/release/1.79.0/source/boost_1_79_0.tar.gz
tar -xvf boost_1_79_0.tar.gz

echo "Building boost"
mkdir -p ${cwd}/boost
cd boost_1_79_0
./bootstrap.sh --prefix=${cwd}/boost/
./b2 --with-iostreams --with-program_options install

cd $cwd
mkdir -p build
cd build

echo "Running CMAKE"
cmake .. -DBOOST_ROOT=${cwd}/boost/

echo "Compiling code"
make -j 8
cd $cwd
