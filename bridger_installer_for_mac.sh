#!/bin/bash

# brewのルートパスを取得
BREW_PATH=`which brew | sed 's/\/bin\/brew//'`

# brewを確認
if [ ! "$BREW_PATH" ] ; then
	echo "brew not installed"
	exit 1
fi

# llvmを確認
if [ ! -d "$BREW_PATH/opt/llvm" ] ; then
	echo "llvm not installed"
	exit 1
fi

# boostを確認
if [ ! -d "$BREW_PATH/include/boost" ] ; then
	echo "boost not installed"
	exit 1
fi

# パスを追加
PATH="$BREW_PATH/opt/llvm/bin:$PATH"

# ファイル書き換え
sed -e 's/return aNode(edge)/return this->aNode(edge)/' -e 's/return bNode(edge)/return this->bNode(edge)/' src/lemon/bits/base_extender.h >src/lemon/bits/base_extender.h~
mv src/lemon/bits/base_extender.h~ src/lemon/bits/base_extender.h

sed -e 's/#include <boost\/graph\/adjacency_list.hpp>/#include <boost\/graph\/vector_as_graph.hpp>\
#include <boost\/graph\/adjacency_list.hpp>/' src/path_search.cpp >src/path_search.cpp~
mv src/path_search.cpp~ src/path_search.cpp

sed -e 's/#include <bits\/typesizes.h>//' -e 's/#include <bits\/types.h>//' src/common.h >src/common.h~
mv src/common.h~ src/common.h

sed -e 's/-R\$boost_ldpath/-rpath,\$boost_ldpath/g' configure >configure~
mv configure~ configure

sed -e 's/\$FindBin::Bin/\$FindBin::RealBin/' -e 's/ directory\// \$output_directory\//' -e 's/print "\\n### Splicing Graphs Reconstruction ###\\n\\n";/print "\\n### Splicing Graphs Reconstruction ###\\n\\n";\
    mkdir "RawGraphs" or die "Error, cannot make \$output_directory\/RawGraphs!\\n";/' Bridger.pl >Bridger.pl~
mv Bridger.pl~ Bridger.pl

chmod 755 Bridger.pl configure

# コンパイル
export CC=clang
export CXX=clang++
./configure --prefix=$PWD --enable-intel64 --with-boost=$BREW_PATH && make AM_MAKEFLAGS='CC=clang CXX=clang++'
mkdir bin && cd bin && ln -s ../Bridger.pl
