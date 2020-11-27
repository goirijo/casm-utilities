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

#git submodule update --init --recursive

echo "Appending plugins to Makefiles"
include_sockets
concat_configures

autoreconf -ifv

