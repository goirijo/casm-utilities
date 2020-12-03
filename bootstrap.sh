include_sockets()
{
    target=plugins/Makesocket.am
    echo "" > ${target}
    for f in $(find -L plugins -mindepth 2 -maxdepth 2 -name Makesocket.am); do
        echo "include ${f}" >> ${target}
    done
}

concat_configures()
{
    target=./configure.ac
    cp _configure.ac ${target}
    for f in $(find -L plugins -mindepth 2 -maxdepth 2 -name configure.ac); do
        cat $f >> ${target}
    done

    echo "AC_OUTPUT" >> ${target}
}

git submodule update --init --recursive

#Automatically generate Makefullcasm.am to track changes from submodules/CASMcode
cd submodules/CASMcode
bash bootstrap.sh > /dev/null
cd ../../
echo "Updating Makefile in casmblob..."
cp submodules/CASMcode/src/casm/Makemodule.am lib/casmblob/Makefullcasm.am
sed -i 's/lib_LTLIBRARIES/noinst_LTLIBRARIES/g' lib/casmblob/Makefullcasm.am
sed -i 's/src\/casm/submodules\/CASMcode\/src\/casm/g' lib/casmblob/Makefullcasm.am
sed -i 's/include\/casm/submodules\/CASMcode\/include\/casm/g' lib/casmblob/Makefullcasm.am
sed -i 's/libcasm/libcasmblob/g' lib/casmblob/Makefullcasm.am


echo "Appending plugins to Makefiles"
include_sockets
concat_configures

autoreconf -ifv

