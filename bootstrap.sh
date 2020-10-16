include_sockets()
{
    target=plugins/Makesocket.am
    echo "" > ${target}
    for f in $(find plugins -mindepth 2 -maxdepth 2 -name Makesocket.am); do
        echo "include ${f}" >> ${target}
    done
}

echo "Appending plugins to Makefiles"
include_sockets

autoreconf -ifv

