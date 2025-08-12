rm -r build
mkdir -p build
pushd build
  cmake ../
  make
popd

